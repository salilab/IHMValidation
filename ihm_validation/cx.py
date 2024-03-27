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
import ihm
import numpy as np
from pathlib import Path
from collections import defaultdict
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

pd.options.mode.chained_assignment = None
NA = 'Not available'


# asym_id, seq_id, atom_id
def get_hierarchy_from_atoms(atoms):
    def infinite_defaultdict(): return defaultdict(infinite_defaultdict)
    root = infinite_defaultdict()

    for a in atoms:
        root[a.asym_unit.id][a.seq_id][a.atom_id] = a

    return root


# asym_id, seq_id, atom_id
def get_hierarchy_from_model(model):
    def infinite_defaultdict(): return defaultdict(infinite_defaultdict)
    root = infinite_defaultdict()

    for a in model.get_atoms():
        root[a.asym_unit.id][a.seq_id][a.atom_id] = a

    for r in model.representation:
        if r.granularity == 'by-residue':
            for i in range(r.asym_unit.seq_id_range[0],
                           r.asym_unit.seq_id_range[1] + 1):
                root[r.asym_unit.asym.id][i]['CA'] = None

    for s in model.get_spheres():
        # Consider only by-residue spheres
        seq_ids = list(set(s.seq_id_range))
        if len(seq_ids) != 1:
            continue

        seq_id = seq_ids[0]

        if root[s.asym_unit.id][seq_id]['CA'] is None:
            root[s.asym_unit.id][seq_id]['CA'] = s

    return root


class CxValidation(GetInputInformation):
    ID = None
    driver = None

    def __init__(self, mmcif_file, cache):
        super().__init__(mmcif_file)
        self.ID = str(self.get_id())
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
                r1n = exl.residue1.comp.id
                r2n = exl.residue2.comp.id

                # Select atoms
                if xl.granularity == 'by-atom':
                    a1n = xl.atom1
                    a2n = xl.atom2

                elif xl.granularity == 'by-residue':
                    # Crosslink applied to the specific residue
                    # represented by the alpha carbon atom
                    a1n = 'CA'
                    a2n = 'CA'
                else:
                    logging.debug(Exception('Unsupported xl granularity'))
                    continue

                r1 = xl.asym1.residue(exl.residue1.seq_id)
                r2 = xl.asym2.residue(exl.residue2.seq_id)

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
        chid = row['chain1']
        rid = row['resnum1']
        an = row['name1']

        a1 = model[chid][rid][an]

        if type(a1) not in allowed_particle_types:
            a1 = None
        if a1 is None:
            logging.debug(f'Atom {chid} {rid} {an} is empty')

        chid = row['chain2']
        rid = row['resnum2']
        an = row['name2']

        a2 = model[chid][rid][an]

        if type(a2) not in allowed_particle_types:
            a2 = None
        if a2 is None:
            logging.debug(f'Atom {chid} {rid} {an} is empty')

        if a1 is None or a2 is None:
            d = None
        else:
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
                cmp = np.isclose(
                    ed, threshold,
                    rtol=row['psi'] + row['sigma1'] + row['sigma2']
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
                    title = f'Satisfaction rates in model group {gimg}'
                    p.title.text = title

                    p.title.text_font_size = "12pt"
                    p.xaxis.axis_label_text_font_size = "14pt"
                    p.yaxis.axis_label_text_font_size = "14pt"
                    p.xaxis.major_label_text_font_size = "14pt"
                    p.yaxis.major_label_text_font_size = "14pt"

                    col = gridplot(
                        [p], ncols=1, toolbar_location='right',
                        # sizing_mode='scale_width'
                    )
                    tab = Panel(child=col, title=f'Model group {gimg}')
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

                p.output_backend = "svg"

                title = f"Model group {gimg}; {lt}: {rt}, {d:.1f} Å"

                p.title.text_font_size = "12pt"
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"

                p.ray(
                    x=d, y=0,
                    line_color='black', angle=np.pi / 2,
                    line_width=2
                )
                p.xaxis.axis_label = 'Euclidean distance, Å'
                p.yaxis.axis_label = 'Count'
                p.title.text = title
                plots.append(p)

            col = gridplot(
                plots, ncols=1, toolbar_location='right',
                # sizing_mode='scale_width'
            )
            tab = Panel(
                child=col, title=f'Model group {gimg}',)
            tabs_.append(tab)

        tabs = Tabs(tabs=tabs_,
                    # sizing_mode='scale_width'
                    )

        title = 'cx_distograms'
        return self.save_plots(tabs, title, imgDirname)

    def save_plots(self, plot, title, imgDirname='.'):
        stem = f'{self.ID}_{title}'

        imgpath = Path(
            imgDirname,
            f'{stem}.html')
        save(
            plot, imgpath,
            resources=CDN,
            title='Satisfaction rates per ensemble',
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
                   webdriver=self.driver)

        svgs = [Path(x).name for x in svgs]

        return (imgpath, imgpath_json, svgs)

    @staticmethod
    def request_pride(url: str) -> dict:
        ''' pull data from PRIDE using crosslinking PDB-DEV API '''
        result = None
        r = requests.get(url)

        if not r.ok:
            # Wait until cold request completes and go to DB cache
            logging.info('Retrying pulling data from PRIDE')
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
        url = f"https://www.ebi.ac.uk/pride/ws/archive/crosslinking/pdbdev/projects/{pid}/sequences"
        result = self.request_pride(url)['data']
        return result

    def get_residue_pairs_pride(self, pid: str, page_size: int = 99) -> dict:
        '''get sequences from PRIDE entry'''
        url = f"https://www.ebi.ac.uk/pride/ws/archive/crosslinking/pdbdev/projects/{pid}/residue-pairs/based-on-reported-psm-level/passing"
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
        cache_fn = Path(self.cache, f'{code}.pickle')
        data = None

        # Check if we already requested the data
        if Path(cache_fn).is_file():
            logging.info(f'Found {code} in cache! {cache_fn}')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)
        elif not Path(cache_fn).is_file():
            ms_seqs = self.get_sequences_pride(code)
            ms_res_pairs = self.get_residue_pairs_pride(code)

            data = {
                'pride_id': code,
                'sequences': ms_seqs,
                'residue_pairs': ms_res_pairs
            }

            with open(cache_fn, 'wb') as f:
                pickle.dump(data, f)

        return data

    def get_pride_ids(self) -> list:
        '''
        get a list of all PRIDE ids from entry
        '''
        ids = []
        for i, d in enumerate(self.system.orphan_datasets):
            if isinstance(d, ihm.dataset.MassSpecDataset) or isinstance(d, ihm.dataset.CXMSDataset):
                if isinstance(d.location, ihm.location.PRIDELocation):
                    try:
                        pid = d.location.access_code
                        ids.append(pid)
                    except AttributeError:
                        pass

        return ids

    def validate_pride_data(self, data: dict) -> dict:

        out = None

        # Unpack pride data
        pid = data['pride_id']
        ms_res_pairs = data['residue_pairs']
        ms_seqs_ = data['sequences']

        if len(ms_seqs_) == 0 or len(ms_res_pairs) == 0:
            return out

        # Get sequences from mmcif entry
        mmcif_seqs = {}
        for e in self.system.entities:
            if e.is_polymeric:
                seq = ''.join([x.code_canonical for x in e.sequence])
                seq = pyhmmer.easel.TextSequence(sequence=seq, name=e._id.encode('utf-8'), description=e.description.encode('utf-8')).digitize(pyhmmer.easel.Alphabet.amino())
                mmcif_seqs[e._id] = seq

        # Get sequences from crosslinking-MS data
        ms_seqs = [
            pyhmmer.easel.TextSequence(
                sequence=x['sequence'],
                name=x['id'].encode('utf-8'),
                source=x['file'].encode('utf-8')
            ).digitize(pyhmmer.easel.Alphabet.amino()) for x in ms_seqs_
        ]

        matches_seqs = {}
        matches_seqs_ids = {}
        matches_residue_pairs = {}
        for k, v in mmcif_seqs.items():
            best_hit = list(pyhmmer.hmmer.phmmer(v, ms_seqs))[0][0]
            matches_seqs[k] = best_hit
            matches_seqs_ids[k] = best_hit.best_domain.hit.name.decode("utf-8")

        matched_ms_seqs = list(matches_seqs_ids.values())
        matched_mmcif_entities = list(matches_seqs_ids.keys())

        ms_rps_total_ = len(ms_res_pairs)
        ms_rps_mmcif_entities = 0

        sxls = []
        for r in ms_res_pairs:
            if r['prot1'] in matched_ms_seqs and r['prot2'] in matched_ms_seqs:
                eid1 = r['prot1']
                rid1 = r['pos1']
                eid2 = r['prot2']
                rid2 = r['pos2']
                sxl = tuple(sorted(((eid1, rid1), (eid2, rid2))))
                sxls.append(sxl)

        sxls = set(sxls)

        ms_rps_mmcif_entities = len(sxls)

        exls = []

        for restr_ in self.system.restraints:
            # We are interested only in Chemical crosslinks
            if type(restr_) != ihm.restraint.CrossLinkRestraint:
                continue

            # Iterate over all crosslinks in the dataset
            for xl in restr_.cross_links:
                eid1 = xl.experimental_cross_link.residue1.entity._id
                rid1 = xl.experimental_cross_link.residue1.seq_id
                eid2 = xl.experimental_cross_link.residue1.entity._id
                rid2 = xl.experimental_cross_link.residue2.seq_id
                exl = tuple(sorted(((eid1, rid1), (eid2, rid2))))
                exls.append(exl)

        exls = list(set(exls))


        mmcif_rps_ms_entities = 0

        emxls = []
        for (eid1, rid1), (eid2, rid2) in exls:
            if eid1 in matched_mmcif_entities and eid2 in matched_mmcif_entities:
                mmcif_rps_ms_entities += 1
                eid1_ = matches_seqs_ids[eid1]
                eid2_ = matches_seqs_ids[eid2]
                exl = tuple(sorted(((eid1_, rid1), (eid2_, rid2))))
                emxls.append(exl)

        rps_both = len(set(emxls) & set(sxls))
        mmcif_rps = len(exls)
        ms_rps = len(ms_res_pairs)

        out = {
            'pride_id': pid,
            'entities_ms': len(ms_seqs),
            'entities': len(mmcif_seqs),
            'matches': [
                {
                    'entity': k,
                    'entity_ms': v.best_domain.hit.name.decode("utf-8"),
                    'e-value': v.best_domain.c_evalue,
                    'exact_match': v.best_domain.alignment.target_sequence == v.best_domain.alignment.hmm_sequence.upper()
                } for k, v in matches_seqs.items()
            ],
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

        return out

    def validate_all_pride_data(self) -> list:
        '''perform data quality validation for all crosslinking-MS datasets'''

        codes = self.get_pride_ids()
        outs = []
        for code in codes:
            data = self.get_pride_data(code)
            out = self.validate_pride_data(data)
            if out is not None:
                outs.append(out)

        return outs
