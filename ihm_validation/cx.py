###################################
# Script :
# 1) Contains class for XL-MS validation
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

from mmcif_io import GetInputInformation
from utility import get_hierarchy_from_model, NA
import pandas as pd
import logging
import ihm
import numpy as np
from pathlib import Path
from bokeh.plotting import save
from bokeh.layouts import gridplot
from bokeh.models import Range1d
from bokeh.models.widgets import Panel, Tabs
from bokeh.resources import CDN
import iqplot
import json
from bokeh.embed import json_item
from bokeh.io import export_svgs
import requests
import pickle
import pyhmmer
import time
import utility
import xml.etree.ElementTree as ET

pd.options.mode.chained_assignment = None


class CxValidation(GetInputInformation):
    ID = None
    driver = None

    def __init__(self, mmcif_file, cache):
        super().__init__(mmcif_file)
        self.cache = cache
        self.nos = self.get_number_of_models()
        self.dataset = self.get_dataset_comp()
        # Only atomic structures are supported so far
        # self.struct = prody.parseMMCIF(mmcif_file, header=False)

        self.entities = set([x.details for x in self.system.asym_units])
        self.chains = set([x.id for x in self.system.asym_units])
        self.get_cx_data()

    def select_atom_by_asym_id_seq_id_atom_id(atoms, asym_id, seq_id, atom_id):
        result = None
        for a in atoms:
            if a.asym_unit.id == asym_id:
                if a.seq_id == seq_id:
                    if a.atom_id == atom_id:
                        result = a
        return result

    def get_models(self):

        # parse all models
        models = {}
        gim = 0
        for istg, stg in enumerate(self.system.state_groups):
            for ist, st in enumerate(stg):
                for img, mg in enumerate(st):
                    for im, m in enumerate(mg):
                        gim += 1
                        m_ = get_hierarchy_from_model(m)
                        models[gim] = m_
        return models

    def get_raw_restraints(self):
        """Get all restraints"""
        allr = []

        rid = -1
        # Iterate over all restraints datasets
        for restr_ in self.system.restraints:
            # We are interested only in Chemical crosslinks
            if type(restr_) != ihm.restraint.CrossLinkRestraint:
                continue

            # Iterate over all crosslinks in the dataset
            for xl in restr_.cross_links:
                # we have same ID for measurments in multiple models
                rid += 1

                # get corresponding experimental crosslink
                exl = xl.experimental_cross_link

                # Extract residue names from atoms
                try:
                    r1n = exl.residue1.comp.id
                    r2n = exl.residue2.comp.id
                except IndexError as e:
                    logging.error('Missing residue')
                    logging.error(e)
                    continue

                # Select atoms
                if xl.granularity == 'by-atom':
                    a1n = xl.atom1
                    a2n = xl.atom2

                elif xl.granularity == 'by-residue':
                    # Crosslink applied to the specific residue
                    # represented by the alpha carbon atom
                    a1n = 'CA'
                    a2n = 'CA'
                elif xl.granularity == 'by-feature':
                    a1n = 'coarse-grained'
                    a2n = 'coarse-grained'
                else:
                    logging.debug(Exception('Unsupported xl granularity'))
                    continue

                try:
                    r1 = xl.asym1.residue(exl.residue1.seq_id)
                    r2 = xl.asym2.residue(exl.residue2.seq_id)
                except IndexError as e:
                    logging.error('Missing residue')
                    logging.error(e)
                    continue

                intra_chain = False
                if xl.asym1.id == xl.asym2.id:
                    intra_chain = True

                intra_entity = False
                if r1.asym.entity.description == r2.asym.entity.description:
                    intra_entity = True

                r_ = {
                    'chemistry': restr_.linker.auth_name,
                    'restraint_id': int(xl._id),
                    'group_id': int(exl._id),
                    'chain1': xl.asym1.id,
                    'resnum1': r1.seq_id,
                    'resnum1_auth': r1.auth_seq_id,
                    'resname1': r1n,
                    'name1': a1n,
                    'chain2': xl.asym2.id,
                    'resnum2': r2.seq_id,
                    'resnum2_auth': r2.auth_seq_id,
                    'resname2': r2n,
                    'name2': a2n,
                    'distance_limit': xl.distance.distance,
                    'distance_lower_limit': xl.distance.distance_lower_limit,
                    'distance_upper_limit': xl.distance.distance_upper_limit,
                    'restraint_type': xl.distance.restraint_type,
                    'psi': xl.psi,
                    'sigma1': xl.sigma1,
                    'sigma2': xl.sigma2,
                    'group_restraint_all': xl.restrain_all,
                    'entity_name1': r1.asym.entity.description,
                    'entity_name2': r2.asym.entity.description,
                    'intra_chain': intra_chain,
                    'intra_entity': intra_entity,
                    # New property to distinguish restraint types
                    'restraint_enum': None,
                    # New property to select restraint group/threshold types
                    'restraint_rtd': None,
                    # Geom properties
                    'state_group': None,
                    'state': None,
                    'model_group': None,
                    'model_number': None,
                    'distance_euclidean': None,
                }

                allr.append(r_)

        # Convert data to pandas dataframe
        allr = pd.DataFrame(allr)

        return allr

    def get_rtdtype(self, row):
        rt = row['restraint_type']

        lt = None

        if row['intra_entity']:
            lt = 'Self-links'
        else:
            lt = 'Heteromeric links'

        if rt == 'upper bound':
            d = row['distance_upper_limit']
        elif rt == 'lower bound':
            d = row['distance_lower_limit']
        elif rt == 'harmonic':
            d = row['distance_limit']
        else:
            raise ValueError('Wrong restraint type')

        return (rt, d, lt)

    def get_ertype(self, row):
        r1n = row['resname1']
        a1n = row['name1']
        r2n = row['resname2']
        a2n = row['name2']

        # Sort names of residues and atoms
        (r1n_, a1n_), (r2n_, a2n_) = sorted(
            [(r1n, a1n), (r2n, a2n)]
        )

        # Construct the extended restraint type
        rtype_ = (
            row['chemistry'],
            r1n_, a1n_, r2n_, a2n_,
            row['restraint_type'],
            # row['distance_lower_limit'],
            row['distance_limit'],
            # row['distance_upper_limit'],
         )

        return rtype_

    def get_ertypes(self):
        ertypes = {}  # enumerate restraints types
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self.get_ertype(row)

            if rtype_ not in ertypes:
                ertypes[rtype_] = len(ertypes)

        return ertypes

    def get_rtdtypes(self):
        ertypes = {}  # enumerate restraints types
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self.get_rtdtype(row)

            if rtype_ not in ertypes:
                ertypes[rtype_] = len(ertypes)

        return ertypes

    def assign_ertypes(self):
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self.get_ertype(row)
            ertype_ = self.ertypes[rtype_]
            self.raw_restraints.at[index, 'restraint_enum'] = ertype_

    def assign_rtdtypes(self):
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self.get_rtdtype(row)
            ertype_ = self.rtdtypes[rtype_]
            self.raw_restraints.at[index, 'restraint_rtd'] = ertype_

    def get_measured_restraints(self):
        restraints = []

        gistg = 0
        gist = 0
        gimg = 0
        gim = 0
        for istg, stg in enumerate(self.system.state_groups):
            gistg += 1
            for ist, st in enumerate(stg):
                gist += 1
                for img, mg in enumerate(st):
                    gimg += 1
                    for im, m in enumerate(mg):
                        gim += 1
                        m_ = self.models[gim]
                        logging.info(f'Assessing crosslinking-MS for MODEL {gim}')
                        for index, row in self.raw_restraints.iterrows():
                            d = self.measure_restraint(m_, row)

                            ndata = {
                                # Store as much information as we can
                                'distance_euclidean': d,
                                'model_number': gim,
                                'model_group': gimg,
                                'state': gist,
                                'state_group': gistg,
                            }

                            nrow = row.to_dict()
                            nrow.update(ndata)
                            restraints.append(nrow)

        restraints = pd.DataFrame(restraints)

        return restraints

    def measure_restraint(self, model, row):
        allowed_particle_types = (ihm.model.Atom, ihm.model.Sphere)
        # Check that we have all necessary atoms
        rid = row['restraint_id']
        gid = row['group_id']
        chid = row['chain1']
        rid = row['resnum1']
        an = row['name1']

        a1 = model[chid][rid][an]

        if not isinstance(a1, allowed_particle_types):
            a1 = None
        if a1 is None:
            logging.warning(f'Restraint {rid}: Atom {chid} {rid} {an} is empty')

        chid = row['chain2']
        rid = row['resnum2']
        an = row['name2']

        a2 = model[chid][rid][an]

        if not isinstance(a2, allowed_particle_types):
            a2 = None
        if a2 is None:
            logging.warning(f'Restraint {rid}: Atom {chid} {rid} {an} is empty')

        if a1 is None or a2 is None:
            d = None
        elif row['name1'] == 'coarse-grained' or row['name2'] == 'coarse-grained':

            if row['name1'] == row['name2'] == 'coarse-grained':
                # Calculate distance between spheres
                a1_ = np.array([a1.x, a1.y, a1.z])
                r1_ = a1.radius
                a2_ = np.array([a2.x, a2.y, a2.z])
                r2_ = a2.radius
                # In case spheres overlap
                d = max(0, np.linalg.norm(a2_ - a1_) - (r1_ + r2_))

            else:
                logging.warning(r'Incompatible crosslinking-MS granularities')
                d = None
        else:
            # Assume atomic distances
            # Calculate distance
            a1_ = np.array([a1.x, a1.y, a1.z])
            a2_ = np.array([a2.x, a2.y, a2.z])
            d = np.linalg.norm(a2_ - a1_)

        return d

    def quality_check(self, data: pd.DataFrame) -> None:
        """Check consistency of crosslink restraints"""
        gids = list(set(data['group_id']))

        for gid in gids:
            data_ = data[data['group_id'] == gid]
            self.check_conditional_flag(data_)

    def check_conditional_flag(self, data: pd.DataFrame) -> None:
        """Check consistency of conditional flags in a restraint group"""
        gid = list(set(data['group_id']))[0]
        conditional_flags = list(set(data['group_restraint_all']))

        if len(conditional_flags) != 1:
            raise ValueError(
                f'Mixed conditional flags in crosslink restraint group {gid}'
            )

    def get_number_of_restraints(self) -> int:
        return len(self.raw_restraints)

    def get_number_of_restraint_groups(self) -> int:
        nrg = len(set(self.raw_restraints['group_id']))
        return nrg

    def get_cx_data(self) -> (pd.DataFrame, pd.DataFrame):
        """Extract crosslinking-MS data from mmcif file"""

        output = (None, None)
        raw_restraints = self.get_raw_restraints()

        if len(raw_restraints) > 0:
            self.quality_check(raw_restraints)
            self.raw_restraints = raw_restraints

            self.ertypes = self.get_ertypes()
            self.assign_ertypes()
            self.ertypes_df = self.get_ertypes_df()

            self.rtdtypes = self.get_rtdtypes()
            self.assign_rtdtypes()

            self.models = self.get_models()
            measured_restraints = self.get_measured_restraints()

            if len(measured_restraints) > 0:

                # Drop missing restraints
                count_missing = measured_restraints[
                    'distance_euclidean'].isna()
                logging.debug(
                    f'Dropped {count_missing} crosslink restraints with '
                    'empty distances'
                )
                measured_restraints.dropna(
                    subset=['distance_euclidean'], inplace=True)

                self.measured_restraints = measured_restraints

                output = (self.ertypes_df, self.measured_restraints)

        return output

    def get_ertypes_df(self):
        # Exctract subset of data as restraint types
        fields_ = [
            'Linker',
            'Residue 1', 'Atom 1',
            'Residue 2', 'Atom 2',
            'Restraint type',
            # 'Lower_limit, Å',
            'Distance, Å',
            # 'Upper_limit, Å'
        ]

        ertypes = pd.DataFrame(
            self.ertypes.keys(),
            index=self.ertypes.values(),
            columns=fields_)

        ertypes['Count'] = None
        for index, row in ertypes.iterrows():
            n = len(
                self.raw_restraints[
                    self.raw_restraints['restraint_enum'] == index]
            )
            ertypes.at[index, 'Count'] = n

        return ertypes

    def get_ertypes_df_html(self):

        return self.ertypes_df.to_dict('tight')

    @staticmethod
    def format_pct_count(pct, count):
        if pct is None:
            pct_ = NA
            pctr_ = NA
        elif isinstance(pct, float):
            pct_ = f'{pct:.2f}'
            pctr_ = f'{100.0 - pct:.2f}'

        else:
            raise Exception('Wrong percentage of violated restraints')

        stats = {
            'Satisfied': pct_,
            'Violated': pctr_,
            'Count': count
            }
        return stats

    def get_best_distance_per_restraint(self, data):
        # Verify, that there is only one type of restraints in the group
        best = None

        rtype = data.iloc[0]['restraint_type']
        dists = data['distance_euclidean'].to_numpy()
        threshold = data.iloc[0]['distance_limit']

        if rtype == 'upper bound':
            best = min(dists)
        elif rtype == 'lower bound':
            diff = (dists - threshold) >= 0

            if diff.any():
                best = min(dists[diff])
            else:
                best = max(dists)
        elif rtype == 'harmonic':
            best = dists[np.argmin(np.abs(dists - threshold))]

        return best

    def get_best_distances_per_model_group(self):
        rtd_groups = {}

        rids = list(set(self.measured_restraints['restraint_id']))

        gistg = 0
        gist = 0
        gimg = 0
        for istg, stg in enumerate(self.system.state_groups):
            gistg += 1
            for ist, st in enumerate(stg):
                gist += 1
                for img, mg in enumerate(st):
                    gimg += 1
                    rtd_groups[gimg] = {}
                    for rid in rids:
                        data_ = self.measured_restraints[
                            (self.measured_restraints['restraint_id'] == rid) &
                            (self.measured_restraints['model_group'] == gimg)]

                        if len(data_) == 0:
                            continue

                        rtdtype = self.get_rtdtype(data_.iloc[0])
                        d = self.get_best_distance_per_restraint(data_)

                        if rtdtype not in rtd_groups[gimg]:
                            rtd_groups[gimg][rtdtype] = []

                        rtd_groups[gimg][rtdtype].append(d)

                    rtd_groups[gimg] = dict(
                        sorted(
                            rtd_groups[gimg].items(),
                            key=lambda x: f'{x[0][2]}_{x[0][0]}_{x[0][1]}'
                        )
                    )

        return rtd_groups

    def is_restraint_group_satisfied(self, data):
        # Verify, that there is only one type of restraints in the group
        satisfied_restraints = np.zeros(len(data), dtype=bool)
        conditional_flag_all = list(set(data['group_restraint_all']))[0]

        for i, (index, row) in enumerate(data.iterrows()):
            rtype = row['restraint_type']
            ed = row['distance_euclidean']
            threshold = row['distance_limit']

            if rtype == 'upper bound':
                cmp = (ed - threshold) <= 0
            elif rtype == 'lower bound':
                cmp = (ed - threshold) >= 0
            # Check with Ben
            elif rtype == 'harmonic':
                atol = 1e-08

                tol_keys = ['sigma1', 'sigma2']

                for k in tol_keys:
                    if row[k] is not None:
                        atol += row[k]

                cmp = np.isclose(
                    ed, threshold,
                    atol=atol
                )

            satisfied_restraints[i] = cmp

        if conditional_flag_all:
            satisfied = satisfied_restraints.all()
        else:
            satisfied = satisfied_restraints.any()
        return satisfied

    def get_restraint_group_chain_type(self, data):
        rg_type = None
        if data['intra_chain'].all():
            rg_type = 'Intramolecular'
        elif (~data['intra_chain']).all():
            rg_type = 'Intermolecular'
        elif data['intra_chain'].any() and (~data['intra_chain']).any():
            rg_type = 'Ambiguous'

        return rg_type

    def get_restraint_group_entity_type(self, data):
        rg_type = None
        if data['intra_entity'].all():
            rg_type = 'Self-links'
        elif (~data['intra_entity']).all():
            rg_type = 'Heteromeric links'
        elif data['intra_entity'].any() and (~data['intra_entity']).any():
            rg_type = 'Ambiguous entity'

        return rg_type

    def process_restraint_groups(self, data, mode='entity'):
        rgs = set(data['group_id'])

        stats = {'All': {'satisfied': 0, 'total': 0}}

        for rgs_ in rgs:
            good_ = 0
            rg_chain_type = None
            rg_entity_type = None
            data_ = data[data['group_id'] == rgs_]
            rg_chain_type = self.get_restraint_group_chain_type(data_)
            rg_entity_type = self.get_restraint_group_entity_type(data_)

            assert rg_chain_type is not None
            assert rg_entity_type is not None

            rg_type = f'{rg_entity_type}/{rg_chain_type}'

            # Check per model
            models = list(set(data_['model_number']))
            for model in models:
                data__ = data_[data_['model_number'] == model]

                if self.is_restraint_group_satisfied(data__):
                    good__ = 1
                else:
                    good__ = 0

                good_ = good_ or good__

            if rg_type not in stats:
                stats[rg_type] = {'satisfied': 0, 'total': 0}

            stats[rg_type]['satisfied'] += good_
            stats[rg_type]['total'] += 1

            stats['All']['satisfied'] += good_
            stats['All']['total'] += 1

        return stats

    def get_stats_per_model_group(self):
        stats = {}

        gistg = 0
        gist = 0
        gimg = 0
        # gim = 0
        for istg, stg in enumerate(self.system.state_groups):
            gistg += 1
            stats[gistg] = {}
            for ist, st in enumerate(stg):
                gist += 1
                stats[gistg][gist] = {}
                for img, mg in enumerate(st):
                    gimg += 1
                    stats[gistg][gist][gimg] = {
                        'cx_stats': None,
                        'ens_stats': None}

                    cx_stats = {}

                    data_ = self.measured_restraints[
                        self.measured_restraints['model_group'] == gimg]

                    stats_ = self.process_restraint_groups(data_)

                    for k, v in stats_.items():
                        pct = None
                        count = v['total']

                        if count > 0:
                            pct = v['satisfied'] / v['total'] * 100.0

                        cx_stats[k] = self.format_pct_count(pct, count)

                    stats[gistg][gist][gimg]['cx_stats'] = cx_stats

                    ens_stats = {}

                    found_ensemble = False
                    for i, e in enumerate(self.system.ensembles, 1):
                        if e.model_group == mg:
                            found_ensemble = True
                            break

                    if found_ensemble:
                        try:
                            ens_stats['ensemble_id'] = i
                            ens_stats['num_models_deposited'] = \
                                e.num_models_deposited
                            ens_stats['num_models'] = e.num_models
                        except AttributeError:
                            logging.error(f'Ens: {i} | Missing model_group')
                        except TypeError:
                            ens_stats['ensemble_id'] = i
                            # Wait for fix in the python-ihm
                            ens_stats['num_models_deposited'] = len(mg)
                            ens_stats['num_models'] = e.num_models
                    else:
                        ens_stats['ensemble_id'] = NA
                        ens_stats['num_models_deposited'] = len(mg)
                        ens_stats['num_models'] = len(mg)

                    stats[gistg][gist][gimg]['ens_stats'] = ens_stats

        return stats


    def get_per_model_satifaction_rates(self) -> list:
        out = []
        gistg = 0
        gist = 0
        gimg = 0
        gim = 0
        tabs_ = []
        for istg, stg in enumerate(self.system.state_groups):
            gistg += 1
            for ist, st in enumerate(stg):
                gist += 1
                for img, mg in enumerate(st):
                    gimg += 1

                    for im, m in enumerate(mg):
                        gim += 1
                        data_ = self.measured_restraints[
                            self.measured_restraints['model_number'] == gim]

                        stats_ = self.process_restraint_groups(data_)

                        for k, v in stats_.items():
                            pct = None
                            count = v['total']

                            if count > 0:
                                pct = v['satisfied'] / v['total'] * 100.0
                                # out_stats[k].append(pct)
                                r = {
                                    'model_number': gim,
                                    'Category': k,
                                    'Satisfaction': pct
                                }

                                if k == 'All':
                                    out.append(pct)

        return(out)

    def plot_satisfaction_per_ensemble(self, imgDirname='.'):
        def scatter_plot(stats):
            xmin, xmax = -2, 102

            jitter = None
            if len(stats['model_number'].unique()) > 1:
                jitter = 'jitter'

            p = iqplot.stripbox(
                data=stats, q='Satisfaction', cats='Category', spread=jitter,
                frame_width=350,
                frame_height=max(100, len(set(stats['Category'])) * 50),
                marker_kwargs=dict(alpha=0.5, size=7)
            )

            p.output_backend = "svg"

            p.x_range = Range1d(xmin, xmax)
            p.xaxis.axis_label = 'Satisfaction rate, %'

            return p

        gistg = 0
        gist = 0
        gimg = 0
        gim = 0
        tabs_ = []
        for istg, stg in enumerate(self.system.state_groups):
            gistg += 1
            for ist, st in enumerate(stg):
                gist += 1
                for img, mg in enumerate(st):
                    gimg += 1

                    out_stats_ = []

                    for im, m in enumerate(mg):
                        gim += 1
                        data_ = self.measured_restraints[
                            self.measured_restraints['model_number'] == gim]

                        stats_ = self.process_restraint_groups(data_)

                        for k, v in stats_.items():
                            pct = None
                            count = v['total']

                            if count > 0:
                                pct = v['satisfied'] / v['total'] * 100.0
                                # out_stats[k].append(pct)
                                r = {
                                    'model_number': gim,
                                    'Category': k,
                                    'Satisfaction': pct
                                }

                                out_stats_.append(r)

                    if len(out_stats_) == 0:
                        continue

                    out_stats = pd.DataFrame(out_stats_)

                    p = scatter_plot(out_stats)
                    title = f'Satisfaction rates in Model Group {gimg}'
                    p.title.text = title

                    p.title.text_font_size = "12pt"
                    p.xaxis.axis_label_text_font_size = "14pt"
                    p.yaxis.axis_label_text_font_size = "14pt"
                    p.xaxis.major_label_text_font_size = "14pt"
                    p.yaxis.major_label_text_font_size = "14pt"
                    p.yaxis.major_label_text_align = 'right'
                    p.yaxis.group_text_align = 'right'
                    p.yaxis.subgroup_text_align = 'right'
                    p.min_border_bottom = 75

                    col = gridplot(
                        [p], ncols=1, toolbar_location='right',
                        # sizing_mode='scale_width'
                    )
                    tab = Panel(child=col, title=f'Model Group {gimg}')
                    tabs_.append(tab)

        tabs = Tabs(tabs=tabs_)

        title = 'cx_ensemble_satisfaction'
        return self.save_plots(tabs, title, imgDirname)

    def plot_distograms_per_model_group(self, imgDirname='.'):
        """plot all restraints in the dataset"""

        data = self.get_best_distances_per_model_group()
        tabs_ = []
        for gimg in data.keys():
            plots = []
            if len(data[gimg]) == 0:
                continue

            for (rt, d, lt), dists in data[gimg].items():
                data_ = pd.DataFrame(dists, columns=['Crosslinks'])

                xmax = 1
                bins = [0, 1]
                if len(dists) > 0:
                    xmax = int(np.ceil(max(dists)))
                    bins = np.linspace(0, xmax, xmax + 1)
                elif len(dists) > 10:
                    bins = 'freedman-diaconis'

                p = iqplot.histogram(
                    data=data_, q='Crosslinks', density=False,
                    bins=bins,
                    style="step_filled",
                    frame_width=500, frame_height=100,
                    # sizing_mode='scale_width',
                )
                p.yaxis.ticker.desired_num_ticks = 3

                p.output_backend = "svg"

                title = f"Model Group {gimg}; {lt}: {rt}, {d:.1f} Å"

                p.title.text_font_size = "12pt"
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_align = 'right'

                p.ray(
                    x=d, y=0,
                    line_color='black', angle=np.pi / 2,
                    line_width=2
                )
                p.xaxis.axis_label = 'Euclidean distance, Å'
                p.yaxis.axis_label = 'Count'
                p.title.text = title
                p.min_border_bottom = 75
                plots.append(p)

            col = gridplot(
                plots, ncols=1, toolbar_location='right',
                # sizing_mode='scale_width'
            )
            tab = Panel(
                child=col, title=f'Model Group {gimg}',)
            tabs_.append(tab)

        tabs = Tabs(tabs=tabs_,
                    # sizing_mode='scale_width'
                    )

        title = 'cx_distograms'
        return self.save_plots(tabs, title, imgDirname)

    def save_plots(self, plot, title, imgDirname='.'):
        stem = f'{self.ID_f}_{title}'

        imgpath = Path(
            imgDirname,
            f'{stem}.html')
        save(
            plot, imgpath,
            resources=CDN,
            title=title,
        )

        imgpath_json = Path(
            imgDirname,
            f'{stem}.json')

        with open(imgpath_json, 'w') as f:
            json.dump(json_item(plot, title), f)

        imgpath_svg = Path(
            imgDirname,
            f'{stem}.svg')

        svgs = export_svgs(plot, filename=imgpath_svg,
                   webdriver=self.driver, timeout=15)

        svgs = [Path(x).name for x in svgs]

        return (imgpath, imgpath_json, svgs)

    @staticmethod
    def request_pride(url: str) -> dict:
        ''' pull data from PRIDE using crosslinking PDB-IHM API '''
        result = None
        r = requests.get(url)

        if not r.ok:
            # Wait until cold request completes and go to DB cache
            logging.info(f'Retrying pulling {url} from PRIDE')
            time.sleep(60)
            r = requests.get(url)

        if r.ok:
            try:
                result = r.json()
            except JSONDecodeError:
                pass

        return result

    def get_sequences_pride(self, pid: str) -> dict:
        '''get sequences from PRIDE entry'''
        result = None
        url = f"https://www.ebi.ac.uk/pride/ws/archive/crosslinking/v2/pdbdev/projects/{pid}/sequences"
        result = self.request_pride(url)
        return result

    def get_residue_pairs_pride(self, pid: str, page_size: int = 99) -> dict:
        '''get sequences from PRIDE entry'''
        url = f"https://www.ebi.ac.uk/pride/ws/archive/crosslinking/v2/pdbdev/projects/{pid}/residue-pairs/based-on-reported-psm-level/passing"
        page = 1
        url_ = f"{url}?page={page}&page_size={page_size}"
        result = self.request_pride(url_)

        rps = []
        if result is not None and 'page' in result:
            session = requests.Session()

            max_page = int(result['page']["total_pages"])
            page_size = int(result['page']["page_size"])

            rps = []
            for i in range(1, max_page + 1):
                url_ = f"{url}?page={i}&page_size={page_size}"
                rps_ = session.get(url_).json()['data']
                rps.extend(rps_)

        return rps

    def get_pride_data(self, code):
        '''
        get data from PRIDE
        '''
        cache_fn = Path(self.cache, f'{code}.pkl')
        data = None

        # Check if we already requested the data
        if Path(cache_fn).is_file():
            logging.info(f'Found {code} in cache! {cache_fn}')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)
        elif not Path(cache_fn).is_file():
            ms_seqs = self.get_sequences_pride(code)
            ms_res_pairs = self.get_residue_pairs_pride(code)

            if ms_seqs is not None and len(ms_res_pairs) > 0:
                data = {
                    'pride_id': code,
                    'sequences': ms_seqs,
                    'residue_pairs': ms_res_pairs
                }

                with open(cache_fn, 'wb') as f:
                    pickle.dump(data, f)

            else:
                logging.info(f'PRIDE data for {code} is incomplete')

        return data

    def get_pride_ids(self) -> list:
        '''
        get a list of all PRIDE ids from entry
        '''
        ids = []
        for i, d in enumerate(self.system.orphan_datasets):
            if isinstance(d, ihm.dataset.CXMSDataset):
                if isinstance(d.location, ihm.location.PRIDELocation) or \
                        isinstance(d.location, ihm.location.ProteomeXchangeLocation):
                    try:
                        pid = d.location.access_code
                        ids.append(pid)
                    except AttributeError:
                        pass
                # Try to automatically convert jPOST ids to PRIDE ids
                if isinstance(d.location, ihm.location.JPOSTLocation):
                    try:
                        pid_ = d.location.access_code
                        r = requests.get(f'https://repository.jpostdb.org/xml/{pid_}.0.xml')
                        xml = ET.fromstring(r.content)
                        pid = xml.find('Project').attrib['pxid']
                        ids.append(pid)
                        logging.info(f'Found PRIDE ID {pid} for JPOST ID {pid_}')
                    # blanket catch because there are too many
                    # potential network exceptions
                    except Exception as e:
                        logging.error(e)
                        pass

        return ids

    def validate_pride_data(self, data: dict) -> tuple :
        """Match sequences, residues pairs and return stats"""
        out = (None, None, None)

        # Unpack pride data
        pid = data['pride_id']
        ms_res_pairs = data['residue_pairs']
        ms_seqs_ = data['sequences']

        if len(ms_seqs_) == 0 or len(ms_res_pairs) == 0:
            return out

        # Get sequences from mmcif entry
        mmcif_seqs = {}
        mmcif_seqs_descriptions = {}

        for e in self.system.entities:
            if e.is_polymeric:
                seq = ''.join([x.code_canonical for x in e.sequence])
                desc = e.description
                seq_ = pyhmmer.easel.TextSequence(
                            sequence=seq,
                            name=e._id.encode('utf-8'),
                            description=desc.encode('utf-8')
                            ).digitize(pyhmmer.easel.Alphabet.amino())

                mmcif_seqs[e._id] = seq_
                mmcif_seqs_descriptions[e._id] = desc

        # Get sequences from crosslinking-MS data
        ms_seqs = [
            pyhmmer.easel.TextSequence(
                sequence=x['sequence'],
                name=x['id'].encode('utf-8'),
                source=x['file'].encode('utf-8')
            ).digitize(pyhmmer.easel.Alphabet.amino()) for x in ms_seqs_
        ]

        # Match sequences using pyHMMER
        # select only 1st best match
        matched_seqs = {}
        matched_seqs_mapping = {}
        matched_seqs_ids = {}

        for k, v in mmcif_seqs.items():
            matches_ = list(pyhmmer.hmmer.phmmer(v, ms_seqs))[0]
            if len(matches_) > 0:
                best_hit = list(pyhmmer.hmmer.phmmer(v, ms_seqs))[0][0]
                matched_seqs[k] = best_hit
                mapping_, _ = self.pyhmmer_alignment_to_map(best_hit)
                matched_seqs_mapping[k] = mapping_
                matched_seqs_ids[k] = best_hit.best_domain.hit.name.decode("utf-8")
            else:
                logging.warning(f"Couldn't match mmCIF entity {k} to any entities in {pid}")
                matched_seqs[k] = None
                matched_seqs_ids[k] = None

        matched_mmcif_entities = list(matched_seqs_ids.keys())
        matched_ms_seqs = list(matched_seqs_ids.values())

        # Filter residue pairs from MS data
        ms_rps_filtered = []
        for r in ms_res_pairs:
            eid1 = r['prot1']
            rid1 = r['pos1']
            eid2 = r['prot2']
            rid2 = r['pos2']
            sxl = tuple(sorted(((eid1, rid1), (eid2, rid2))))
            ms_rps_filtered.append(sxl)

        ms_rps_filtered = set(ms_rps_filtered)

        # Select MS residue pairs from matched sequences
        ms_rps_mmcif_entities = 0
        sxls = []
        for r in ms_rps_filtered:
            (eid1, rid1), (eid2, rid2) = r

            if eid1 in matched_ms_seqs and eid2 in matched_ms_seqs:
                sxls.append(r)

        sxls = set(sxls)
        ms_rps_mmcif_entities = len(sxls)

        # Select residue pairs from the entry
        exls = []

        for restr_ in self.system.restraints:
            # We are interested only in Chemical crosslinks
            if type(restr_) != ihm.restraint.CrossLinkRestraint:
                continue

            # Iterate over all crosslinks in the dataset
            for xl in restr_.cross_links:
                eid1 = xl.experimental_cross_link.residue1.entity._id
                rid1 = xl.experimental_cross_link.residue1.seq_id
                eid2 = xl.experimental_cross_link.residue2.entity._id
                rid2 = xl.experimental_cross_link.residue2.seq_id

                exl = tuple(sorted(((eid1, rid1), (eid2, rid2))))
                exls.append(exl)

        exls = list(set(exls))

        # Find corresponding entry - MS data crosslinks
        mmcif_rps_ms_entities = 0

        rps_mapping = []
        emxls = []
        for rps in exls:
            (eid1, rid1), (eid2, rid2) = rps

            if eid1 in matched_mmcif_entities and eid2 in matched_mmcif_entities:
                mmcif_rps_ms_entities += 1
                eid1_ = matched_seqs_ids[eid1]
                eid2_ = matched_seqs_ids[eid2]

                try:
                    rid1_ = matched_seqs_mapping[eid1][rid1]
                except KeyError:
                    logging.debug(f"Can't map residue {eid1} {rid1} to {eid1_}")
                    continue

                try:
                    rid2_ = matched_seqs_mapping[eid2][rid2]
                except KeyError:
                    logging.debug(f"Can't map residue {eid2} {rid2} to {eid2_}")
                    continue

                if eid1_ is None or eid2_ is None:
                    continue

                exl = tuple(sorted(((eid1_, rid1_), (eid2_, rid2_))))
                if exl in sxls:
                    # Good matching crosslinks
                    rps_mapping.append((rps, exl))
                else:
                    # Residue pairs unique to the entry
                    rps_mapping.append((rps, None))

                emxls.append(exl)

        emxls = set(emxls)

        # Add non-mapped residue pairs from MS data
        for exl in list(sxls):
            if exl not in emxls:
                rps_mapping.append((None, exl))

        # Calculate some stats
        rps_both = len(set(emxls) & set(sxls))
        mmcif_rps = len(exls)
        ms_rps = len(ms_rps_filtered)

        # Prepare output

        out = {
            'pride_id': pid,
            'entities_ms': len(ms_seqs),
            'entities': len(mmcif_seqs),
            'matches': [],
            'stats': {
                'entry': {
                    'total': mmcif_rps,
                    'mapped_entities': mmcif_rps_ms_entities,
                    'mapped_entities_pct': mmcif_rps_ms_entities / mmcif_rps * 100.,
                    'matched': rps_both,
                    'matched_pct': rps_both / mmcif_rps * 100.,
                },
                'ms': {
                    'total': ms_rps,
                    'mapped_entities': ms_rps_mmcif_entities,
                    'mapped_entities_pct': ms_rps_mmcif_entities / ms_rps * 100.,
                    'matched': rps_both,
                    'matched_pct': rps_both / ms_rps * 100.,
                }
            }
        }

        # Add stats about matches
        for k, v in matched_seqs.items():
            if v is not None:
                match_ =  {
                    'entity': k,
                    'entity_desc': mmcif_seqs_descriptions[k],
                    'entity_ms': v.best_domain.hit.name.decode("utf-8"),
                    'e-value': v.best_domain.c_evalue,
                    'exact_match': v.best_domain.alignment.target_sequence == v.best_domain.alignment.hmm_sequence.upper(),
                }

            else:
                match_ =  {
                    'entity': k,
                    'entity_desc': mmcif_seqs_descriptions[k],
                    'entity_ms': utility.NA,
                    'e-value': utility.NA,
                    'exact_match': utility.NA,
                }


            out['matches'].append(match_)

        return (out, matched_seqs, rps_mapping)

    def validate_all_pride_data(self) -> list:
        '''perform data quality validation for all crosslinking-MS datasets'''

        codes = self.get_pride_ids()
        outs = []
        for code in codes:
            data = self.get_pride_data(code)
            if data is not None:
                out, _, __ = self.validate_pride_data(data)
                if out is not None:
                    outs.append(out)

        return outs

    @staticmethod
    def get_pyhmmer_version():
        """return pyhmmer version"""
        return pyhmmer.__version__

    @staticmethod
    def pyhmmer_alignment_to_map(hit) -> (dict, list):
        """Convert HMMER alignment into residue map"""

        # gaps in HMMER text format
        GAPS = set(['-', '.'])

        mapping_raw = []
        mapping_short = {}
        aln = hit.best_domain.alignment
        mmcif_start = aln.hmm_from
        mmcif_seq = aln.hmm_sequence
        db_start = aln.target_from
        db_seq = aln.target_sequence

        ii = mmcif_start - 1
        ij = db_start - 1

        # iterate over alignment
        # residue indices start from 1
        for i, (aai, aaj) in enumerate(zip(mmcif_seq.upper(), db_seq.upper())):
            if aai not in GAPS:
                ii += 1

            if aaj not in GAPS:
                ij += 1

            if len(set([aai, aaj]) & GAPS) == 0:
                mapping_raw.append(((ii, aai), (ij, aaj)))
                mapping_short[ii] = ij

        # return dict to map residue indices and raw mapping data
        return mapping_short, mapping_raw
