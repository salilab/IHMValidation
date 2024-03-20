#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Validation framework
This is a beta-version of a new validation framework
"""
import logging
from typing import Union
from collections import defaultdict

import ihm
import numpy as np
import pandas as pd

class Validator(object):
    '''Base validator class'''
    dataset = None
    restraint = None

    def __init__(self):
        pass

    def load_restraint(self, restraint: ihm.restraint.Restraint):
        '''Parse the restraint data'''
        raise NotImplementedError(
            'Method not implemented for that class')

    def validate_model(self, model: ihm.model.Model) -> dict:
        '''Validate the model against the data'''
        raise NotImplementedError(
            'Model validation is not implemented for that class')

    def validate_ensemble(self, ensemble) -> dict:
        '''Validate the ensemble against the data'''
        raise NotImplementedError(
            'Ensemble validation is not implemented for that class')

class CXMSValidator(Validator):
    def __init__(self):
        super().__init__()
        self.supported_particles = (ihm.model.Atom, ihm.model.Sphere)

    def load_restraint(self, restraint: ihm.restraint.CrossLinkRestraint) -> None:
        """Extract crosslinking-MS data from mmcif file"""

        raw_restraints = self._load_restraint_raw(restraint)

        if len(raw_restraints) > 0:
            self._quality_check(raw_restraints)
            self.raw_restraints = raw_restraints
            self._assign_ertypes()
            self._assign_rtdtypes()

    def _load_restraint_raw(self, restraint: ihm.restraint.CrossLinkRestraint):
        """Get raw restraint data"""
        allr = []

        # Iterate over all crosslinks in the dataset
        for xl in restraint.cross_links:
            exl = xl.experimental_cross_link

            # Get residues
            r1 = xl.asym1.residue(exl.residue1.seq_id)
            r2 = xl.asym2.residue(exl.residue2.seq_id)

            # Get residue names
            r1n = r1.comp.id
            r2n = r2.comp.id

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

            intra_chain = False
            if xl.asym1.id == xl.asym2.id:
                intra_chain = True

            intra_entity = False
            if r1.asym.entity.description == r2.asym.entity.description:
                intra_entity = True

            r_ = {
                'chemistry': restraint.linker.auth_name,
                'restraint_id': xl.id,
                'group_id': xl.group_id,
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

    def _quality_check(self, data: pd.DataFrame) -> None:
        """Check consistency of crosslink restraints"""
        gids = list(set(data['group_id']))

        for gid in gids:
            data_ = data[data['group_id'] == gid]
            self._check_conditional_flag(data_)

    def _check_conditional_flag(self, data: pd.DataFrame) -> None:
        """Check consistency of conditional flags in a restraint group"""
        gid = list(set(data['group_id']))[0]
        conditional_flags = list(set(data['group_restraint_all']))

        if len(conditional_flags) != 1:
            raise ValueError(
                f'Mixed conditional flags in crosslink restraint group {gid}'
            )

    @property
    def number_of_restraints(self) -> int:
        return len(self.raw_restraints)

    @property
    def number_of_restraint_groups(self) -> int:
        nrg = len(set(self.raw_restraints['group_id']))
        return nrg


    def _get_rtdtype(self, row):
        '''Determine crosslink types: Self/Heteromeric; Upper/Lower/Harmonic; Threshold'''
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

    def _get_ertype(self, row):
        '''Get extended restraint type, that includes linker, residue names, atom names'''
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

    @property
    def ertypes(self):
        '''Get all extendend restraint types'''
        ertypes = {}  # enumerate restraints types
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self._get_ertype(row)

            if rtype_ not in ertypes:
                ertypes[rtype_] = len(ertypes)

        return ertypes

    @property
    def rtdtypes(self):
        '''Get all crosslink types'''
        ertypes = {}  # enumerate restraints types
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self._get_rtdtype(row)

            if rtype_ not in ertypes:
                ertypes[rtype_] = len(ertypes)

        return ertypes

    def _assign_ertypes(self):
        '''Assign extended restraint types'''
        ertypes = self.ertypes
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self._get_ertype(row)
            ertype_ = ertypes[rtype_]
            self.raw_restraints.at[index, 'restraint_enum'] = ertype_

    def _assign_rtdtypes(self):
        '''Assign crosslink types'''
        rtdtypes = self.rtdtypes
        for index, row in self.raw_restraints.iterrows():
            rtype_ = self._get_rtdtype(row)
            ertype_ = rtdtypes[rtype_]
            self.raw_restraints.at[index, 'restraint_rtd'] = ertype_

    @property
    def ertypes_df(self):
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

    def _measure_restraints(self, models: list) -> pd.DataFrame:
        restraints = []
        for m in models:
            m_ = get_hierarchy_from_model(m)
            for index, row in self.raw_restraints.iterrows():
                d = self._measure_restraint(m_, row)

                ndata = {
                    # Store as much information as we can
                    'distance_euclidean': d,
                    'model_number': int(m._id),
                }

                nrow = row.to_dict()
                nrow.update(ndata)
                restraints.append(nrow)

        restraints = pd.DataFrame(restraints)


        if len(restraints) > 0:

            # Drop missing restraints
            count_missing = restraints[
                'distance_euclidean'].isna()
            logging.debug(
                f'Dropped {count_missing} crosslink restraints with '
                'empty distances'
            )
            restraints.dropna(
                subset=['distance_euclidean'], inplace=True)

        return restraints

    def _measure_restraint(self, model, row) -> pd.DataFrame:
        # Check that we have all necessary atoms
        chid = row['chain1']
        rid = row['resnum1']
        an = row['name1']

        a1 = model[chid][rid][an]

        if type(a1) not in self.supported_particles:
            a1 = None
        if a1 is None:
            logging.debug(f'Atom {chid} {rid} {an} is empty')

        chid = row['chain2']
        rid = row['resnum2']
        an = row['name2']

        a2 = model[chid][rid][an]

        if type(a2) not in self.supported_particles:
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

    def validate_model(self, model: ihm.model.Model) -> dict:

        data = self._measure_restraints([model])
        out = self._process_satisfaction_data(data)

        return(out)

    def validate_ensemble(self, models: list) -> dict:

        data = self._measure_restraints(models)
        out = self._process_satisfaction_data(data)
        return(out)

    def _process_satisfaction_data(self, data: pd.DataFrame) -> dict:
        out = {}
        stats = self._process_restraint_groups(data)

        for k, v in stats.items():
            pct = None
            count = v['Count']

            if count > 0:
                pct = v['Satisfied'] / v['Count'] * 100.0
                out[k] = {'Satisfaction': pct, 'Count': count}

        return out

    def _process_restraint_groups(self, data: pd.DataFrame, mode='entity') -> dict:
        rgs = set(data['group_id'])

        stats = {'All': {'Satisfied': 0, 'Count': 0}}

        for rgs_ in rgs:
            good_ = 0
            rg_chain_type = None
            rg_entity_type = None
            data_ = data[data['group_id'] == rgs_]
            rg_chain_type = self._get_restraint_group_chain_type(data_)
            rg_entity_type = self._get_restraint_group_entity_type(data_)

            assert rg_chain_type is not None
            assert rg_entity_type is not None

            rg_type = f'{rg_entity_type}/{rg_chain_type}'

            # Check per model
            models = list(set(data_['model_number']))
            for model in models:
                data__ = data_[data_['model_number'] == model]

                if self._is_restraint_group_satisfied(data__):
                    good__ = 1
                else:
                    good__ = 0

                good_ = good_ or good__

            if rg_type not in stats:
                stats[rg_type] = {'Satisfied': 0, 'Count': 0}

            stats[rg_type]['Satisfied'] += good_
            stats[rg_type]['Count'] += 1

            stats['All']['Satisfied'] += good_
            stats['All']['Count'] += 1

        return stats

    def _is_restraint_group_satisfied(self, data: pd.DataFrame) -> bool:
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

    def _get_restraint_group_chain_type(self, data: pd.DataFrame) -> str:
        rg_type = None
        if data['intra_chain'].all():
            rg_type = 'Intramolecular'
        elif (~data['intra_chain']).all():
            rg_type = 'Intermolecular'
        elif data['intra_chain'].any() and (~data['intra_chain']).any():
            rg_type = 'Ambiguous'

        return rg_type

    def _get_restraint_group_entity_type(self, data: pd.DataFrame) -> str:
        rg_type = None
        if data['intra_entity'].all():
            rg_type = 'Self-links'
        elif (~data['intra_entity']).all():
            rg_type = 'Heteromeric links'
        elif data['intra_entity'].any() and (~data['intra_entity']).any():
            rg_type = 'Ambiguous entity'

        return rg_type


class StereoChemistryValidator(Validator):
    def __init__(self):
        super().__init__()
    def validate_model(self, model: ihm.model.Model) -> dict:
        '''Validate the model against the data'''
        if is_model_atomic(model):
            pass
        pass


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

def is_model_mixed(model: ihm.model.Model) -> bool:
    """Check if model is atomic"""
    result = False
    granularities = set([r.granularity for r in model.representation])
    if len(granularities) > 1:
        result = True
    return result

def is_model_atomic(model: ihm.model.Model) -> bool:
    """Check if model is atomic"""
    result = False
    granularities = set([r.granularity for r in model.representation])
    if granularities == set(['by-atom']):
        result = True
    return result

def is_model_cg(model: ihm.model.Model) -> bool:
    """Check if model is atomic"""
    result = False
    granularities = set([r.granularity for r in model.representation])
    if granularities == set(['by-residue', 'by-feature']):
        result = True
    return result

