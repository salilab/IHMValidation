# -*- coding: utf-8 -*-
#
# em.py - 3DEM validation for PDB-IHM
#
# Copyright (C) 2025 Arthur Zalevsky <aozalevsky@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
3DEM validation for PDB-IHM
"""
# system
import logging
import json
import pickle
import io
import re
import requests
import time
import gzip
import subprocess

# handling files
from pathlib import Path
import shutil
import glob
import tempfile
import os

# numerical
import numpy as np
import pandas as pd

# plotting
from bokeh.plotting import save, figure
from bokeh.layouts import gridplot
from bokeh.models import Range1d
from bokeh.models.widgets import Panel, Tabs
from bokeh.embed import json_item
from bokeh.io import export_svgs
import iqplot

# ihmv
from mmcif_io import GetInputInformation
import utility
import ihm
from molprobity import MyModelDumper, AtomSiteVariant, _stub
import images

# EMDB
import va

# Set global pandas options
pd.options.mode.chained_assignment = None

class EMValidation(GetInputInformation):
    ID = None
    driver = None
    VA_TIMEOUT = 900

    def __init__(self, mmcif_file, cache):
        super().__init__(mmcif_file)
        self.cache = cache
        # Only atomic structures are supported so far

    def save_plots(self, plot, title, imageDirName='.'):
        """Save bokeh plots as svg images and json blobs"""
        stem = f'{self.ID_f}_{title}'

        imgpath_json = Path(
            imageDirName,
            f'{stem}.json')

        with open(imgpath_json, 'w') as f:
            json.dump(json_item(plot, title), f)

        imgpath_svg = Path(
            imageDirName,
            f'{stem}.svg')

        svgs = export_svgs(plot, filename=imgpath_svg,
                   webdriver=self.driver, timeout=15)

        svgs = [Path(x).name for x in svgs]

        return ('', imgpath_json, svgs)

    @staticmethod
    def request_emdb(url: str) -> dict:
        '''pull data from EMDB'''
        result = None
        r = requests.get(url)

        if not r.ok:
            # Wait until cold request completes and go to DB cache
            logging.info(f'Retrying pulling {url} from EMDB')
            time.sleep(60)
            r = requests.get(url)

        if r.ok:
            try:
                result = r.json()
            except JSONDecodeError:
                pass

        return result

    def get_emdb_data(self, code):
        '''
        get data from EMDB
        '''
        cache_fn = Path(self.cache, f'{code}.pkl')
        data = None

        # Check if we already requested the data
        if Path(cache_fn).is_file():
            logging.info(f'Found {code} in cache! {cache_fn}')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)
        elif not Path(cache_fn).is_file():
            map_metadata = self.get_emdb_map_metadata(code)
            map_validation = self.get_emdb_map_validation(code)
            map_path = Path(self.cache, f'{code}.map')
            map_path = self.get_emdb_map(code, map_path)
            map_plots = self.get_emdb_map_images(code, map_validation)

            if (map_metadata is not None and
                map_validation is not None and
                Path(map_path).is_file() and
                len(map_plots) > 0):

                data = {
                    'emdb_id': code,
                    'map_metadata': map_metadata,
                    'map_validation': map_validation,
                    'map_path': map_path,
                    'map_plots': map_plots
                }

                with open(cache_fn, 'wb') as f:
                    pickle.dump(data, f)
            else:
                logging.info(f'EMDB data for {code} is incomplete')

        return data

    @staticmethod
    def get_emdb_numerical_code(code) -> str:
        """Extract only numerical part of the EMDB ID"""
        code_ = re.search('\d+', code).group()
        return code_

    def get_emdb_map(self, code, fname) -> str:
        '''download validation data from EMDB'''
        # Extract only numerical part of the ID
        code_ = self.get_emdb_numerical_code(code)
        out = requests.get(f'https://files.wwpdb.org/pub/emdb/structures/EMD-{code_}/map/emd_{code_}.map.gz', stream=True)
        if out.status_code == 200:
            stream_ = io.BytesIO()
            for chunk in out.raw.stream(1024, decode_content=False):
                stream_.write(chunk)
            data_ = gzip.decompress(stream_.getvalue())
            with open(fname, 'wb') as f:
                f.write(data_)

        else:
            logging.error(f'Could not download EMDB map, status code: {out.status_code}')

        return fname

    def get_emdb_map_metadata(self, code) -> dict:
        '''download validation data from EMDB'''
        out = {}
        # Extract only numerical part of the ID
        code_ = self.get_emdb_numerical_code(code)
        r = requests.get(f'https://www.ebi.ac.uk/emdb/api/entry/{code_}')
        if r.status_code == 200:
            out = r.json()
        else:
            logging.error(f"Couldn't get map metadata from EMDB for {code}")
        return out

    def get_emdb_map_images(self, code, map_validation) -> dict:
        plots = {}

        try:
             plots_ = {}
             for p in ['x', 'y', 'z']:
                 name = f"{p}projection"
                 label = p.upper()
                 plots_[label] = utility.encode_img(self.get_map_image(code, name))
             plots['map_proj'] = plots_
        except (KeyError, IndexError, ValueError, TypeError) as e:
            logging.error(e)

        try:
            plots_ = {}
            for p in ['x', 'y', 'z']:
                name = f"{p}central_slice"
                ind = map_validation['central_slice']['indices'][p]
                label = f"{p.upper()} Index: {ind}"
                plots_[label] = utility.encode_img(self.get_map_image(code, name))
            plots['map_central_slice'] = plots_
        except (KeyError, IndexError, ValueError, TypeError) as e:
            logging.error(e)

        try:
             plots_ = {}
             for p in ['x', 'y', 'z']:
                 name = f"{p}largestvariance_slice"
                 ind = map_validation['largest_variance_slice']['indices'][p]
                 label = f"{p.upper()} Index: {ind}"
                 plots_[label] = utility.encode_img(self.get_map_image(code, name))
             plots['map_largest_variance_slice'] = plots_
        except (KeyError, IndexError, ValueError, TypeError) as e:
            logging.error(e)

        try:
            plots_ = {}
            for p in ['x', 'y', 'z']:
                name = f"glow_{p}std"
                label = p.upper()
                plots_[label] = utility.encode_img(self.get_map_image(code, name))
            plots['map_glow_std'] = plots_
        except (KeyError, IndexError, ValueError, TypeError) as e:
            logging.error(e)

        try:
            plots_ = {}
            for p in ['x', 'y', 'z']:
                name = f"scaled_{p}surface"
                label = p.upper()
                plots_[label] = utility.encode_img(self.get_map_image(code, name))
            plots['map_surface'] = plots_
        except (KeyError, IndexError, ValueError, TypeError) as e:
            logging.error(e)

        # Try rawmap, not always available
        if 'rawmap_orthogonal_projection' in map_validation:
            try:
                plots_ = {}
                for p in ['x', 'y', 'z']:
                    name = f"{p}projection"
                    plots_[p.upper()] = utility.encode_img(self.get_map_image(code, name, rawmap=True))
                plots['rawmap_proj'] = plots_
            except (KeyError, IndexError, ValueError, TypeError) as e:
                logging.error(e)

            try:
                plots_ = {}
                for p in ['x', 'y', 'z']:
                    name = f"{p}central_slice"
                    ind = map_validation['rawmap_central_slice']['indices'][p]
                    label = f"{p.upper()} Index: {ind}"
                    plots_[label] = utility.encode_img(self.get_map_image(code, name, rawmap=True))
                plots['rawmap_central_slice'] = plots_
            except (KeyError, IndexError, ValueError, TypeError) as e:
                logging.error(e)

            try:
                plots_ = {}
                for p in ['x', 'y', 'z']:
                    name = f"{p}largestvariance_slice"
                    ind = map_validation['rawmap_largest_variance_slice']['indices'][p]
                    label = f"{p.upper()} Index: {ind}"
                    plots_[label] = utility.encode_img(self.get_map_image(code, name, rawmap=True))
                plots['rawmap_largest_variance_slice'] = plots_
            except (KeyError, IndexError, ValueError, TypeError) as e:
                logging.error(e)

            try:
                plots_ = {}
                for p in ['x', 'y', 'z']:
                    name = f"glow_{p}std"
                    label = p.upper()
                    plots_[label] = utility.encode_img(self.get_map_image(code, name, rawmap=True))
                plots['rawmap_glow_std'] = plots_
            except (KeyError, IndexError, ValueError, TypeError) as e:
                logging.error(e)

            try:
                plots_ = {}
                for p in ['x', 'y', 'z']:
                    name = f"scaled_{p}surface"
                    label = p.upper()
                    plots_[label] = utility.encode_img(self.get_map_image(code, name, rawmap=True))
                plots['rawmap_surface'] = plots_
            except (KeyError, IndexError, ValueError, TypeError) as e:
                logging.error(e)

        if 'masks' in map_validation:
            plots['masks'] = {}
            for k, v in map_validation['masks'].items():
                try:
                    plots_ = {}
                    for p in ['x', 'y', 'z']:
                        name = f"{k}_scaled_{p}maskview"
                        label = p.upper()
                        plots_[label] = utility.encode_img(self.get_map_image(code, name))
                    plots['masks'][k] = plots_
                except (KeyError, IndexError, ValueError, TypeError) as e:
                    logging.error(e)

        return plots

    def get_emdb_map_validation(self, code) -> dict:
        '''download validation data from EMDB'''
        # Extract only numerical part of the ID
        code_ = self.get_emdb_numerical_code(code)
        r = requests.get(f'https://www.ebi.ac.uk/emdb/api/analysis/{code_}?information=all')
        out = {}
        if r.status_code == 200:
            out = r.json()[code_]
        else:
            logging.error(f"Couldn't get data quality from EMDB for {code}")
        return out

    def get_emdb_ids(self) -> list:
        '''
        get a list of all EMDB ids from entry
        '''
        ids = []
        for i, d in enumerate(self.system.orphan_datasets):
            if isinstance(d, ihm.dataset.EMDensityDataset):
                if isinstance(d.location, ihm.location.EMDBLocation):
                    try:
                        pid = d.location.access_code
                        ids.append(pid)
                    except AttributeError:
                        pass
        return ids

    def get_map_image(self, code, name, rawmap=False, extension="jpeg"):
        img = None
        code_ = self.get_emdb_numerical_code(code)
        url = f"https://www.ebi.ac.uk/pdbe/emdb/emdb-entry/emdbva/va-{code_}/va/emd_{code_}.map_{name}.{extension}"
        if rawmap:
            url = f"https://www.ebi.ac.uk/pdbe/emdb/emdb-entry/emdbva/va-{code_}/va/emd_{code_}_rawmap.map_{name}.{extension}"

        r = requests.get(url, stream=True)
        if r.status_code == 200:
            buffer = io.BytesIO()
            for chunk in r.raw.stream(1024, decode_content=False):
                buffer.write(chunk)

            img = buffer.getvalue()

        else:
            logging.error(f"Couldn't download {url}")

        return img

    def validate_emdb_data(self, data: dict, imageDirName='.') -> tuple :
        """Match sequences, residues pairs and return stats"""
        out = (None, None)
        data_stats, data_plots, data_stats_plots = {}, {}, {}
        fit_stats, fit_plots = {}, {}

        # Unpack emdb data
        emdbid = data['emdb_id']
        emdbid_ = self.get_emdb_numerical_code(emdbid)
        map_metadata = data['map_metadata']
        map_validation = data['map_validation']
        resolution = None
        try:
            resolution = float(map_validation["resolution"]["value"])
        except (KeyError, ValueError, IndexError, TypeError) as e:
            logging.error(e)

        recl = None
        try:
            recl = float(map_validation['recommended_contour_level']['recl'])
        except (KeyError, ValueError, IndexError, TypeError) as e:
            logging.error(e)
        data_plots = data['map_plots']

        # Voxel-value plot
        try:
            p = figure(
                title=f"Voxel-value distribution (Mode={map_validation['density_distribution']['mode']})",
                x_axis_label='Voxel value',
                y_axis_label='Number of voxels (log10)',
                plot_height=350
                )
            X = map_validation['density_distribution']['x']
            Y = map_validation['density_distribution']['y']
            p.line(X, Y, color='blue')
            try:
                p.line([recl, recl], [0, max(Y)], line_color='red', line_width=3, legend_label=f'Recommended contour level {recl:.2f}')
            except (KeyError, ValueError, IndexError, TypeError) as e:
                logging.error(e)
            p.yaxis.ticker.desired_num_ticks = 3
            p.output_backend = "svg"

            p.title.text_font_size = "14pt"
            p.xaxis.axis_label_text_font_size = "14pt"
            p.yaxis.axis_label_text_font_size = "14pt"
            p.xaxis.major_label_text_font_size = "14pt"
            p.yaxis.major_label_text_font_size = "14pt"
            p.xaxis.axis_label_text_font_style = 'normal'
            p.yaxis.axis_label_text_font_style = 'normal'

            title = f'{emdbid}_voxel'
            plots_ = self.save_plots(p, title, imageDirName)
            data_stats_plots[title] = plots_
        except (KeyError, ValueError, IndexError, TypeError) as e:
            logging.error(e)

        # Volume estimate
        vol = None
        try:
            vol = map_validation['volume_estimate']['estvolume']
        except (KeyError, ValueError, IndexError, TypeError) as e:
            logging.error(e)
        try:
            X = map_validation['volume_estimate']['level']
            Y = map_validation['volume_estimate']['volume']

            p = figure(
                    title=f"Volume estimate (Estimated volume={vol:.2f} nm³)",
                    x_axis_label = 'Contour level',
                    y_axis_label = 'Volume [nm³]',
                    plot_height=350
            )
            p.yaxis.ticker.desired_num_ticks = 3
            p.output_backend = "svg"

            p.line(X, Y, color='blue')

            try:
                p.line([recl, recl], [0, max(Y)], color='red', legend_label=f'Recommended contour level {recl:.2f}')
            except (KeyError, ValueError, IndexError, TypeError) as e:
                logging.error(e)

            p.line([min(X), max(X)], [vol, vol], color='orange', legend_label=f'Estimated volume {vol:.2f} nm³')

            title = f'{emdbid}_volume'

            p.title.text_font_size = "14pt"
            p.xaxis.axis_label_text_font_size = "14pt"
            p.yaxis.axis_label_text_font_size = "14pt"
            p.xaxis.major_label_text_font_size = "14pt"
            p.yaxis.major_label_text_font_size = "14pt"
            p.xaxis.axis_label_text_font_style = 'normal'
            p.yaxis.axis_label_text_font_style = 'normal'

            plots_ = self.save_plots(p, title, imageDirName)
            data_stats_plots[title] = plots_
        except (KeyError, ValueError, IndexError, TypeError) as e:
            logging.error(e)


        # RAPS
        if 'rotationally_averaged_power_spectrum' in map_validation:
            try:
                p = figure(
                        title=f"Rotationally averaged power spectrum",
                        x_axis_label = 'Spatial frequency [Å⁻¹]',
                        y_axis_label = 'Log (I)',
                        plot_height=350
                    )
                p.yaxis.ticker.desired_num_ticks = 3
                p.output_backend = "svg"

                X = map_validation['rotationally_averaged_power_spectrum']['x']
                Y = map_validation['rotationally_averaged_power_spectrum']['y']
                maxy = np.max(Y)
                miny = np.min(Y)

                if 'rawmap_rotationally_averaged_power_spectrum' in map_validation:
                    X_raw = map_validation['rawmap_rotationally_averaged_power_spectrum']['x']
                    Y_raw = map_validation['rawmap_rotationally_averaged_power_spectrum']['y']
                    maxy = max(np.max(Y), np.max(Y_raw))
                    miny = min(np.min(Y), np.min(Y))
                    p.line(X_raw, Y_raw, color='orange', legend_label='Raw map RAPS')

                p.line(X, Y, color='blue', legend_label='Primary map RAPS')

                try:
                    loc = 1 / resolution
                    p.line([loc, loc], [miny, maxy], color='red', legend_label=f'Reported resolution {resolution:.2f}*')
                except (KeyError, ValueError, IndexError, TypeError) as e:
                    logging.error(e)

                p.title.text_font_size = "14pt"
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.xaxis.axis_label_text_font_style = 'normal'
                p.yaxis.axis_label_text_font_style = 'normal'

                title = f'{emdbid}_raps'
                plots_ = self.save_plots(p, title, imageDirName)
                data_stats_plots[title] = plots_
            except (KeyError, ValueError, IndexError, TypeError) as e:
                logging.error(e)

        # FSC

        if 'fsc' in map_validation:
            try:
                p = figure(
                    title=f"FSC",
                    x_axis_label = 'Spatial frequency [Å⁻¹]',
                    y_axis_label = 'Correlation',
                    y_range = (-0.1, 1.1),
                    plot_height=350
                    )

                p.yaxis.ticker.desired_num_ticks = 3
                p.output_backend = "svg"

                if 'load_fsc' in map_validation:
                    X0 = map_validation['load_fsc']['curves']['fscx']
                    Y0 = map_validation['load_fsc']['curves']['fscy']
                    p.line(X0, Y0, legend_label='Author provided', color='blue')

                X = map_validation['fsc']['curves']['level']
                Y1 = map_validation['fsc']['curves']['fsc']
                p.line(X, Y1, legend_label='Unmasked-calculated FSC', color='orange')

                Y2 = map_validation['fsc']['curves']['0.143']
                p.line(X, Y2, line_dash='dashed', legend_label='0.143', color='green')

                Y2 = map_validation['fsc']['curves']['0.5']
                p.line(X, Y2, line_dash='dashed', legend_label='0.5', color='red')

                Y3 = map_validation['fsc']['curves']['halfbit']
                p.line(X, Y3, line_dash='dashed', legend_label='Half-bit', color='purple')

                p.line([loc, loc], [miny, maxy], color='black', legend_label=f'Reported resolution {resolution:.2f}*')

                p.title.text_font_size = "14pt"
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.xaxis.axis_label_text_font_style = 'normal'
                p.yaxis.axis_label_text_font_style = 'normal'

                title = f'{emdbid}_fsc'
                plots_ = self.save_plots(p, title, imageDirName)
                data_stats_plots[title] = plots_
            except (KeyError, ValueError, IndexError, TypeError) as e:
                logging.error(e)

        # Fit-to-data
        if self.atomic:
            mid, _, __ = utility.get_larget_assembly_model(self.system)
            mid = int(mid)
            varoot = Path(self.cache, f"{self.ID_f}_va")
            vapath = Path(varoot, f"{mid}_{emdbid}")
            os.makedirs(varoot, exist_ok=True)

            if Path(vapath).is_dir():
                logging.info(f"Found VA for {emdbid} in cache: {vapath}")
            else:
                self.run_va(data, mid, vapath)
            try:
                plots_ = {}
                for p in ['x', 'y', 'z']:
                    fname = glob.glob(str(Path(vapath, f'*cif*.map_scaled_{p}surface.jpeg')))[0]
                    b64 = utility.load_img(fname)
                    plots_[p] = b64
                fit_plots['map_model'] = plots_

                plots_ = {}
                for p in ['x', 'y', 'z']:
                    fname = glob.glob(str(Path(vapath, f'*cif*.map_scaled_{p}qscoresurface.jpeg')))[0]
                    b64 = utility.load_img(fname)
                    plots_[p.upper()] = b64

                fit_plots['model'] = mid
                fit_plots['map_model_q'] = plots_
                fit_plots['q_scale'] = images.QSCORE_SCALE

                plots_ = {}
                for p in ['x', 'y', 'z']:
                    fname = glob.glob(str(Path(vapath, f'*cif*.map_scaled_{p}fitsurface.jpeg')))[0]
                    b64 = utility.load_img(fname)
                    plots_[p] = b64
                fit_plots['map_model_inclusion'] = plots_
                fit_plots['ai_scale'] = images.AI_SCALE

                fname = glob.glob(str(Path(vapath, f'*.map_all.json')))[0]
                emdbid_ = self.get_emdb_numerical_code(emdbid)
                with open(fname, 'r') as f:
                    data_ = json.load(f)[f'EMD-{emdbid_}.map']

                fit_stats[mid] = {'ai_score': {}, 'q_score': {}}
                fit_stats[mid]['ai_score']['chains'] = {}
                fit_stats[mid]['ai_score']['chains']['All'] = {
                    'value': data_['atom_inclusion_by_level']['0']['average_ai_model'],
                    'color': data_['atom_inclusion_by_level']['0']['average_ai_color'],
                    'numberOfAtoms': data_['atom_inclusion_by_level']['0']['totalNumberOfAtoms']
                }
                fit_stats[mid]['ai_score']['chains'].update(data_['atom_inclusion_by_level']['0']['chainaiscore'])
                fit_stats[mid]['ai_score']['average'] = data_['atom_inclusion_by_level']['0']['average_ai_model']

                p = figure(
                    title=f"Atom inclusion",
                    x_axis_label = 'Contour level',
                    y_axis_label = 'Fraction of atoms inside the map',
                    y_range = (0., 1.1),
                    plot_height=350,
                )
                p.output_backend = "svg"

                X = data_['atom_inclusion_by_level']['0']['level']
                Y1 = data_['atom_inclusion_by_level']['0']['backbone']
                Y2 = data_['atom_inclusion_by_level']['0']['all_atom']

                ai_recl_backbone = np.round(Y1[np.abs(recl - np.asarray(X)).argmin()] * 100.).astype(int)
                ai_recl_all = np.round(Y2[np.abs(recl - np.asarray(X)).argmin()] * 100.).astype(int)
                fit_stats[mid]['ai_score']['ai_recl_backbone'] = ai_recl_backbone
                fit_stats[mid]['ai_score']['ai_recl_all'] = ai_recl_all

                p.line(X, Y1, color='blue', legend_label='Backbone atoms')
                p.line(X, Y2, color='orange', legend_label='All non-hydrogen atoms')

                p.line([recl, recl], [0, 1.1], color='red', legend_label=f'Recommended contour level {recl:.3f}')

                p.title.text_font_size = "14pt"
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.xaxis.axis_label_text_font_style = 'normal'
                p.yaxis.axis_label_text_font_style = 'normal'

                title = f'{mid}_{emdbid}_ai_plot'
                plots_ = self.save_plots(p, title, imageDirName)
                fit_plots['ai_plot'] = plots_

                fit_stats[mid]['q_score']['chains'] = {}
                fit_stats[mid]['q_score']['chains']['All'] = {
                    'value': data_['qscore']['0']['data']['averageqscore'],
                    'color': data_['qscore']['0']['data']['averageqscore_color'],
                }
                fit_stats[mid]['q_score']['chains'].update(data_['qscore']['0']['data']['chainqscore'])
                fit_stats[mid]['q_score']['average'] = data_['qscore']['0']['data']['averageqscore']

            except (FileNotFoundError, OSError, IndexError, KeyError) as e:
                logging.error(e)
                fit_stats = {}
                fit_plots = {}

        data_stats['specimen_preparation_list'] = map_metadata['structure_determination_list']['structure_determination'][0]['specimen_preparation_list']['specimen_preparation']
        data_stats['reconstruction_method'] = self.get_em_reconstruction_method(map_metadata)
        resolution_estimates, comments = self.get_resolution_estimates(map_metadata, map_validation)
        data_stats['resolution'] = resolution
        data_stats['resolution_estimates'] = resolution_estimates
        if len(resolution_estimates) == 0:
            data_stats['resolution_estimates'] = None
        data_stats['resolution_estimates_comments'] = comments
        data_stats['recommended_level'] = recl
        if 'rawmap_contour_level' in map_validation:
            data_stats['rawmap_recommended_level'] = float(map_validation['rawmap_contour_level']['cl'])
        data_stats['estimated_volume'] = vol
        out_ = {
            'emdbid': emdbid, 'emdbid_': emdbid_,
            'data_stats': data_stats,
            'data_plots': data_plots,
            'data_stats_plots': data_stats_plots,
            "fit_stats": fit_stats, "fit_plots": fit_plots
        }

        return (out_, None, None)

    def validate_all_emdb_data(self, imageDirName='.') -> list:
        '''perform data quality validation for all crosslinking-MS datasets'''

        codes = self.get_emdb_ids()
        outs = []
        for code in codes:
            data = self.get_emdb_data(code)
            if data is not None:
                out, _, __ = self.validate_emdb_data(data, imageDirName=imageDirName)
                if out is not None:
                    outs.append(out)

        return outs

    @staticmethod
    def get_chimera_version() -> str:
        """return chimera version"""
        version_string = subprocess.check_output(['chimera', '--version', '--nogui']).decode()
        version = re.search(' (\d+.\d+) ', version_string).groups()[0]
        return version

    @staticmethod
    def get_chimerax_version() -> str:
        """return chimera version"""
        version_string = subprocess.check_output(['chimerax', '--version', '--nogui']).decode()
        version = re.search(' (\d+.\d+) ', version_string).groups()[0]
        return version

    @staticmethod
    def get_mapq_version() -> str:
        """return mapq version"""
        with tempfile.NamedTemporaryFile('w') as f:
            f.write('from mapq import mapqVersion; print(mapqVersion)')
            f.flush()
            version_string = subprocess.check_output(['chimera', '--nogui', '--script', f.name]).decode()
            version = version_string.strip()
        return version

    @staticmethod
    def get_va_version() -> str:
        """return va version"""
        version = va.__version__
        return version

    def run_va(self, data, mid, outpath):
        """Run EMDB va validation pipeline for map and model"""

        # Select representative model
        mid = int(mid)

        logging.info(f'Running va from scratch for MODEL {mid}')

        map_path = data['map_path']
        map_metadata = data['map_metadata']
        method = self._get_em_reconstrucion_method_va(map_metadata)
        map_validation = data['map_validation']
        resolution = float(map_validation["resolution"]["value"])
        recl = float(map_validation['recommended_contour_level']['recl'])

        fn = self.mmcif_file
        # tempcif
        outfn = str(Path(self.cache, f'{self.stem}_temp.cif'))

        system, encoding = utility.parse_ihm_cif(fn)

        # Reassign label_asym_id as auth_asym_id
        # for MolProbity
        for asym in system.asym_units:
            asym.auth_seq_id_map = 0
            asym._strand_id = asym.id

        system._check_after_write = _stub

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tempdir:
            outfn = str(Path(tempdir, f'{self.stem}_temp.cif'))

            with open(outfn, 'w', encoding=encoding) as f:
                ihm.dumper.write(f, [system], variant=AtomSiteVariant(mid))

            mapfn = str(Path(map_path).name)
            shutil.copy(map_path, str(Path(tempdir, mapfn)))

            call = [
                # 'xvfb-run', '-a', # Use framebuffer on machines with wayland
                'va',
                '-d', tempdir,
                '-met', str(method),
                '-run', 'surface', 'inclusion', 'qscore',
                '-cl', str(recl),
                '-s', str(resolution),
                '-m', Path(mapfn).name,
                '-f', Path(outfn).name
            ]

            try:
                subprocess.check_call(call, cwd=tempdir, timeout=self.VA_TIMEOUT)
            except subprocess.TimeoutExpired:
                logging.error("Could not finish the VA run")

            shutil.copytree(tempdir, outpath)

    @staticmethod
    def get_level_from_primary_contour(data: dict) -> float:
        """Get recommended level for the map"""
        level = None
        for contour in data['map']['contour_list']['contour']:
            if contour['primary']:
                level = contour['level']

        return level

    @staticmethod
    def get_resolution_estimates(map_data, val_data):
        estimates = {}
        comments = {'general': '', 'fsc': ''}

        resenum = {
            'FSC 0.143 CUT-OFF': '0.143',
            'FSC 0.5 CUT-OFF': '0.5',
            'FSC 1/2 BIT CUT-OFF': 'Half-bit',
            '0.143': '0.143',
            '0.5': '0.5',
            'halfbit': 'Half-bit'
        }

        rresenum = {
            '0.143': 'FSC 0.143 CUT-OFF',
            '0.5': 'FSC 0.5 CUT-OFF',
            'Half-bit': 'FSC 1/2 BIT CUT-OFF',
        }


        try:
            estimates_ = {}
            rec = map_data['structure_determination_list']['structure_determination'][0]['image_processing'][0]['final_reconstruction']
            res = rec['resolution']
            byauthor = res['res_type'] == 'BY AUTHOR'
            if byauthor:
                resmethod = rec['resolution_method']
                if resmethod == 'OTHER':
                    resenum[resmethod] = 'Other'
                if resmethod in resenum:
                    resmethod_ = resenum[resmethod]
                    for k, v in resenum.items():
                        estimates_[v] = '-'
                    estimates_[resmethod_] = f'{float(res["valueOf_"]):.2f}'

                    estimates['Reported by author'] = estimates_
        except KeyError:
            comments['general'] += 'Author provided resolution is not available. '

        if 'load_fsc' in val_data:
            estimates_ = {}
            for k, v in resenum.items():
                estimates_[v] = '-'
            for k, v in val_data['load_fsc']['intersections'].items():
                if k in resenum:
                    estimates_[resenum[k]] = f'{np.round(1 / v["x"], 2):.2f}'

            estimates['Author-provided FSC curve'] = estimates_
        else:
            comments['general'] += 'Author-provided FSC curve is not available. '

        if 'fsc' in val_data:
            estimates_ = {}
            for k, v in resenum.items():
                estimates_[v] = '-'
            for k, v in val_data['fsc']['intersections'].items():
                if k in resenum:
                    estimates_[resenum[k]] = f'{np.round(1 / v["x"], 2):.2f}'
                float(res["valueOf_"])
            estimates['Unmasked-calculated*'] = estimates_
            comments['fsc'] = '*Resolution estimate based on FSC curve calculated by comparison of deposited half-maps. '

            for k, v in estimates_.items():
                ares = estimates['Reported by author'][k]
                if v != '-' and ares != '-':
                    if abs(float(ares) - float(v)) > (float(ares) / 10.):
                        comments['fsc'] += f'The value from deposited half-maps intersecting {rresenum[k]} {v} differs from the reported value {ares} by more than 10%. '

        return estimates, comments

    @staticmethod
    def get_em_reconstruction_method(map_metadata):
        """Get EM reconsruction method from the conrolled vocabulary"""
        m = utility.NA
        methods = {
            'crystallography': 'CRYSTALLOGRAPHY',
            'singleparticle': 'SINGLE PARTICLE',
            'tomography': 'TOMOGRAPHY',
            'subtomogramaveraging': 'SUBTOMOGRAM AVERAGING',
            'helical': 'HELICAL'
        }
        m_ = map_metadata['structure_determination_list']['structure_determination'][0]['method'].lower()
        try:
            m = methods[m_]
        except KeyError as e:
            logging.error(f'Unknown em reconstuction method {e}')

        return m

    @staticmethod
    def _get_em_reconstrucion_method_va(map_metadata):
        """Get EM reconstruction methods for the VA run"""
        m = map_metadata['structure_determination_list']['structure_determination'][0]['method'].lower()
        return m
