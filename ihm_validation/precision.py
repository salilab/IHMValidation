###################################
# Script :
# 1) Contains class for XL-MS validation
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

from mmcif_io import GetInputInformation
import pandas as pd
import logging
from pathlib import Path
import pickle
import time
import utility
# import pymol
import sys

pd.options.mode.chained_assignment = None

p = Path(__file__).absolute().parent.parent.parent
sys.path.append(str(p))
p = Path(p, 'prism', 'src')
sys.path.append(str(p))
p = Path(p, 'prism', 'src')
sys.path.append(str(p))
# Path to pymol
sys.path.append('/usr/lib/python3/dist-packages')

import prism.src.ihm_parser
import prism.src.main
from pymol import cmd

class PRISM(GetInputInformation):
    timeout = 120
    data = {}

    def __init__(self, mmcif_file, cache='.', nocache=False, timeout=120):
        super().__init__(mmcif_file, cache=cache, nocache=nocache)
        self.timeout = timeout
        # Only atomic structures are supported so far
        # self.struct = prody.parseMMCIF(mmcif_file, header=False)

    def get_data(self):
        cache_fn = Path(self.cache, f'{self.stem}.prism.pkl')
        data = {}

        # Check if we already requested the data
        if Path(cache_fn).is_file() and not self.nocache:
            logging.info(f'Found {self.stem} in cache: {cache_fn}')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)

        elif not Path(cache_fn).is_file() or self.nocache:
            logging.info("PrISM analysis is being calculated...")
            data = self._get_data()

            with open(cache_fn, 'wb') as f:
                pickle.dump(data, f)

        self.data = data

        return data

    def _get_data(self):
        data = {}
        for ens in self.system.ensembles:
            mg = ens.model_group

            if mg is None:
                logging.error('ModelGroup is missing for Ensemble {ens._id}. Skipping.')
                continue

            try:
                e_total = ens.num_models
                e_deposited = ens.num_models_deposited
            except (AttributeError, TypeError) as e:
                logging.error('Missing Ensemble {ens._id} attributes')
                logging.error(e)
                continue

            try:
                assert e_deposited == len(mg)
            except AssertionError as e:
                logging.warning('Number of deposited models and size of a model group are different for Ensemble {ens._id}. Using model group size')
                e_deposited = len(mg)

            # Only collect information for multimodel model groups
            if e_deposited < 2:
                continue

            coords, radius, mass, ps_names = prism.src.ihm_parser.get_all_attributes_ihm(mg)
            try:
                # PrISM excecution time is unpredictable
                with utility.timeout(self.timeout):
                    df = prism.src.main.run_prism_ihm( coords, mass, radius, ps_names, classes=3)
            except TimeoutError as e:
                logging.error(f'PrISM calculation timed out for ensemble {ens._id}')
                continue
            except ValueError as e:
                logging.error(f'PrISM calculation failed out for ensemble {ens._id}')
                logging.error(e)
                continue

            data[ens._id] = {
                'total': e_total,
                'deposited': e_deposited,
                'prism': df
            }

        self.data = data

        return data

    def get_plots(self, imgDirname='.'):
        """Render PrISM images with PyMOL"""
        data = self.data
        plots = {}

        def init_colors(cmd):
            # Original PrISM colors

            # cmd.set_color('low_1',  [1.00, 0.00, 0.00])
            # cmd.set_color('low_2',  [1.00, 0.50, 0.50])
            # cmd.set_color('low_3',  [1.00, 0.75, 0.75])
            # cmd.set_color('mid_1',  [1.00, 1.00, 1.00])
            # cmd.set_color('high_3',  [0.56, 0.93, 0.56])
            # cmd.set_color('high_2', [0.20, 0.80, 0.20])
            # cmd.set_color('high_1', [0.00, 0.39, 0.00])

            # https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=7

            cmd.set_color('low_1',  [0.70, 0.09, 0.17])
            cmd.set_color('low_2',  [0.94, 0.54, 0.38])
            cmd.set_color('low_3',  [0.99, 0.86, 0.78])
            cmd.set_color('mid_1',  [0.97, 0.97, 0.97])
            cmd.set_color('high_3', [0.82, 0.89, 0.94])
            cmd.set_color('high_2', [0.40, 0.66, 0.81])
            cmd.set_color('high_1', [0.13, 0.40, 0.68])

        for k, v in data.items():
            cmd.reinitialize()
            init_colors(cmd)
            fn_stem = f'{self.ID_f}_prism_{k}'
            plots_ = []

            for index, row in v['prism'].iterrows():
                x = row['x']
                y = row['y']
                z = row['z']
                r = row['r']
                c = f'{row["Type"]}_{row["Class"]}'
                a = cmd.pseudoatom('', pos=[x, y, z], color=c, vdw=r)

            cmd.show('spheres')
            cmd.bg_color('white')

            cmd.set('opaque_background', 1)
            cmd.orient()
            cmd.zoom(buffer=10)
            cmd.clip('slab', 1000)
            cmd.set('ray_shadow', 0)
            cmd.rotate('z', 180)

            fn = f'{fn_stem}_p1.png'
            cmd.png(str(Path(imgDirname, fn)), ray=False)
            plots_.append(fn)
            cmd.rotate('x', 90)
            fn = f'{fn_stem}_p2.png'
            cmd.png(str(Path(imgDirname, fn)), ray=False)
            plots_.append(fn)
            cmd.rotate('y', 90)
            fn = f'{fn_stem}_p3.png'
            cmd.png(str(Path(imgDirname, fn)), ray=False)
            plots_.append(fn)
            fn = f'{fn_stem}.pse'
            cmd.save(Path(imgDirname, fn))
            plots[k] = {
                'total': v['total'],
                'deposited': v['deposited'],
                'plots': plots_}

        return plots


    @property
    def pymol_version(self):
        return cmd.get_version()[0]

