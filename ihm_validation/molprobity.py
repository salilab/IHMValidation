#!/usr/bin/env python
###################################
# Script :
# 1) Contains class to
# generate molprobity assessments
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import logging
import pickle
import os
from pathlib import Path
import subprocess
from subprocess import run, check_call, CalledProcessError
from mmcif_io import GetInputInformation, MAX_NUM_MODELS
import ihm, ihm.reader, ihm.dumper
import collections
import pandas as pd
import csv
import re
import string
import utility
from typing import Literal
import argparse
import itertools

class GetMolprobityInformation(GetInputInformation):
    _tempfiles = []
    data = {}
    convert = Path(__file__).with_name('molprobity_convert.py')

    def __init__(self, mmcif_file, cache=".", nocache=False):
        super().__init__(mmcif_file, cache, nocache)
        self.verify_molprobity_installation()

        self._tempcif = str(Path(self.cache, f'{self.stem}_temp.cif'))
        # if Path(self._tempcif).is_file():
        #     os.remove(self._tempcif)
        # self.rewrite_mmcif(self._tempcif)
        self._tempfiles.append(self._tempcif)

    def verify_molprobity_installation(self):
        """ Stub for a validation function """
        pass

    @property
    def molprobity_version(self, tool: str = 'molprobity.clashscore') -> str:
        """
        Get MolProbity version.
        We assume that all tools belong to the same release.
        """
        version = None
        try:
            # Try to get "internal version" 4.x.x
            version = self.get_internal_version()
        except OSError:
            # Fallback to commit-based version
            version = subprocess.check_output(
                [tool, '--version'],
                text=True,
                stderr=subprocess.STDOUT).strip()

        return version

    def get_internal_version(self, tool: str = 'molprobity.clashscore') -> str:
        """
        Get internal molprobity version.
        We assume that all tools belong to the same release.
        """
        version = None

        mp_tool_path = subprocess.check_output(
            ['which', tool],
            text=True,
            stderr=subprocess.STDOUT).strip()

        mp_core_path = Path(
            Path(mp_tool_path).parent,
            '../../molprobity/lib/core.php'
        )

        if mp_core_path.is_file():
            with open(str(mp_core_path), 'r') as f:
                raw = f.readlines()

            for line in raw:
               q = re.search('^define\\("MP_VERSION", "(?P<version>.*)"\\);', line)

               if q:
                   version = q.group('version')
                   break

            return version

        else:
            raise OSError('Molprobity core.php module is missing')

    def run_molprobity(self, fname) -> dict|None:
        """ Run MolProbity"""
        mp_stem = Path(fname).with_suffix('.mp.tmp')
        data = None
        # self._tempfiles.append(f_name)

        # Molprobity writes output
        # to prefix.out and prefix.pkl
        try:
            check_call(
                ["molprobity.molprobity",
             "disable_uc_volume_vs_n_atoms_check=True",
             "coot=False",
             "probe=False",
             f"prefix={mp_stem}",
             "pickle=True",
             fname,
             ]
            )
        except CalledProcessError as e:
            logging.error("Couldn't complete MolProbity call")
            logging.error(e)


        else:
            fname_out = f"{mp_stem}.out"
            fname_pkl = f"{mp_stem}.pkl"
            fname_pkl_convert = f"{mp_stem}_convert.pkl"

            run(["mmtbx.python", self.convert, '-i', fname_pkl, '-o', fname_pkl_convert])

            with open(fname_pkl_convert, 'rb') as f:
                data = pickle.load(f)

            os.remove(fname_pkl)
            os.remove(fname_pkl_convert)
            os.remove(fname_out)

        return data

    def get_mp_data(self):

        cache_fn = Path(self.cache, self.stem + '.mp.pkl')

        if cache_fn.is_file() and not self.nocache:
            logging.info(f'Found MolProbity data {cache_fn} in cache')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)

        else:

            logging.info(f'Running MolProbity from scratch')
            fn = self.mmcif_file
            outfn = self._tempcif

            system, encoding = utility.parse_ihm_cif(fn)

            # Reassign label_asym_id as auth_asym_id
            # for MolProbity
            for asym in system.asym_units:
                asym.auth_seq_id_map = 0
                asym._strand_id = asym.id

            system._check_after_write = _stub

            data = {}

            for group, model in system._all_models():

                mid = int(model._id)

                with open(outfn, 'w', encoding=encoding) as f:
                    ihm.dumper.write(f, [system], variant=AtomSiteVariant(mid))
                data_ = self.run_molprobity(outfn)
                if data_ is not None:
                    data[mid] = data_

            os.remove(outfn)

            with open(cache_fn, 'wb') as f:
                pickle.dump(data, f)

        if len(data) == 0:
            logging.warning('Empty MolProbity data')

        self.data = data

        return data

    def summarize_bonds(self):
        data = self.data

        total = 0
        outliers = 0

        duplicates = {}
        for k, v in data.items():
            total += v['bonds']['total']
            outliers += v['bonds']['outliers']

            for r in v['bonds']['outliers_list']:
                a1, a2, observed, ideal, score = r
                key = (a1, a2)
                if key not in duplicates:
                    duplicates[key] = []
                duplicates[key].append((k, r))

        header = (
            'Chain', 'Res', 'Type', 'Atoms',
            '|Z|', 'Observed (Å)', 'Ideal (Å)',
            'Model ID (Worst)', 'Models (Total)'
        )

        duplicates_stats = []

        for k, v in duplicates.items():
            # v is a list of (model_id, data) records
            # data[0:2] - atoms
            assert len(set([x[1][3] for x in v])) == 1
            v_ = sorted(v, key=lambda x: abs(x[1][-1]))

            worst_model_, data_ = v_[-1]

            a1 = data_[0]
            a2 = data_[1]

            chid = a1[0]
            resid = a1[1]
            resname = a1[2]
            name1 = a1[3]
            name2 = a2[3]

            observed_ = data_[2]
            ideal_ = data_[3]
            score_ = data_[4]

            total_ = len(v_)

            r_ = (chid, resid, resname, f"{name1}-{name2}", score_, observed_, ideal_, worst_model_, total_)
            duplicates_stats.append(r_)

        duplicates_stats = sorted(duplicates_stats, key=lambda x: x[4], reverse=True)
        duplicates_stats.insert(0, header)

        return(total, outliers, duplicates_stats)

    def summarize_angles(self):
        data = self.data

        total = 0
        outliers = 0

        duplicates = {}
        for k, v in data.items():
            total += v['angles']['total']
            outliers += v['angles']['outliers']

            for r in v['angles']['outliers_list']:
                a1, a2, a3, observed, ideal, score = r
                key = (a1, a2, a3)
                if key not in duplicates:
                    duplicates[key] = []
                duplicates[key].append((k, r))

        header = (
            'Chain', 'Res', 'Type', 'Atoms',
            '|Z|', 'Observed (Å)', 'Ideal (Å)',
            'Model ID (Worst)', 'Models (Total)'
        )

        duplicates_stats = []

        for k, v in duplicates.items():
            # v is a list of (model_id, data) records
            # data[0:3] - atoms
            try:
                assert len(set([x[1][4] for x in v])) == 1
            except AssertionError as e:
                logging.error('Mixed angle defitions')
                logging.error(e)
                logging.error(v)
                continue

            v_ = sorted(v, key=lambda x: abs(x[1][-1]))
            worst_model_, data_ = v_[-1]

            a1 = data_[0]
            a2 = data_[1]
            a3 = data_[2]

            chid = a1[0]
            resid = a1[1]
            resname = a1[2]
            name1 = a1[3]
            name2 = a2[3]
            name3 = a3[3]

            observed_ = data_[3]
            ideal_ = data_[4]
            score_ = data_[5]

            total_ = len(v_)

            r_ = (
                chid, resid, resname, f"{name1}-{name2}-{name3}",
                score_, observed_, ideal_, worst_model_, total_
            )

            duplicates_stats.append(r_)

        duplicates_stats = sorted(duplicates_stats, key=lambda x: x[4], reverse=True)
        duplicates_stats.insert(0, header)

        return(total, outliers, duplicates_stats)

    def summarize_clashscores(self):
        data = self.data

        header = ('Model ID', 'Clash score', 'Number of clashes')
        stats = []

        for k, v in data.items():
            score_ = v['clash']['clashscore']
            clashes_ = len(v['clash']['clashes_list'])

            r_ = (k, score_, clashes_)

            stats.append(r_)

        stats.insert(0, header)

        return(stats)

    def summarize_clashes(self):
        data = self.data

        total = None
        outliers = 0

        duplicates = {}
        for k, v in data.items():
            outliers += len(v['clash']['clashes_list'])

            for r in v['clash']['clashes_list']:
                a1, a2, score = r
                key = (a1, a2)
                if key not in duplicates:
                    duplicates[key] = []
                duplicates[key].append((k, r))

        header = (
                'Atom 1', 'Atom 2',
                'Clash(Å)',
                'Model ID (Worst)', 'Models (Total)'
            )
        duplicates_stats = []

        for k, v in duplicates.items():
            # v is a list of (model_id, data) records
            # data[0:3] - atoms
            v_ = sorted(v, key=lambda x: abs(x[1][-1]))
            worst_model_, data_ = v_[-1]

            a1 = data_[0]
            a2 = data_[1]

            score_ = data_[2]

            total_ = len(v_)

            r_ = (
                ':'.join(a1), ':'.join(a2), score_,
                worst_model_, total_
            )

            duplicates_stats.append(r_)

        duplicates_stats = sorted(duplicates_stats, key=lambda x: x[2], reverse=True)
        duplicates_stats.insert(0, header)

        return(total, outliers, duplicates_stats)

    def summarize_rama_rota_stats(self, mode: Literal["rama", "rota"]):
        data = self.data

        header = ('Model ID', 'Analysed', 'Favored', 'Allowed', 'Outliers')
        stats = []

        for k, v in data.items():
            total_ = v[mode]['total']
            favored_ = v[mode]['favored']
            allowed_ = v[mode]['allowed']
            outliers_ = v[mode]['outliers']

            r_ = (k, total_, favored_, allowed_, outliers_)
            stats.append(r_)

        stats.insert(0, header)

        return(stats)

    def summarize_rama_rota(self, mode: Literal["rama", "rota"]):
        data = self.data

        duplicates = {}
        for k, v in data.items():
            for r in v[mode]['outliers_list']:
                chid, resid, resname, score_ = r
                key = (chid, resid, resname)
                if key not in duplicates:
                    duplicates[key] = []
                duplicates[key].append((k, r))

        header = (
            'Chain', 'Res', 'Type',
            'Models (Total)'
        )

        duplicates_stats = []

        for k, v in duplicates.items():
            model_id_, data_ = v[0]

            chid = data_[0]
            resid = data_[1]
            resname = data_[2]
            total_ = len(v)

            r_ = (chid, resid, resname, total_)
            duplicates_stats.append(r_)

        duplicates_stats = sorted(duplicates_stats, key=lambda x: (-x[3], x[0], int(x[1])))
        duplicates_stats.insert(0, header)

        return duplicates_stats

    def get_mq_plot_data(self):
        data = self.data

        stats = None

        if len(data) > 0:

            stats = {}

            for k, v in data.items():
                clashscore_ = v['clash']['clashscore']
                rama_ = v['rama']['outliers']
                rota_ = v['rota']['outliers']
                r_ = {
                    'Clashscore': clashscore_,
                    'Ramachandran outliers': rama_,
                    'Sidechain outliers': rota_
                }
                stats[k] = r_

        return stats

    def get_summary_table_stats(self):
        data = self.data

        def get_mp_formatted_range(data, k1, k2, format=".2f"):
            data_ = [x[k1][k2] for x in data.values()]
            r_ = utility.format_range(data_, format=format)

            return r_

        stats = None

        if len(data) > 0:

            stats = [
                f'Clashscore: {get_mp_formatted_range(data, "clash", "clashscore", format=".2f")}',
                f'Ramachandran outliers: {get_mp_formatted_range(data, "rama", "outliers", format="d")}',
                f'Sidechain outliers: {get_mp_formatted_range(data, "rota", "outliers", format="d")}',
            ]

        return stats

    def summarize_mp_data(self):
        mp_data = None

        if len(self.data) > 0:
            mp_data = {}

            total, outliers, outliers_list = self.summarize_bonds()
            mp_data['bonds'] = {}
            mp_data['bonds']['total'] = total
            mp_data['bonds']['outliers'] = outliers
            mp_data['bonds']['outliers_list'] = outliers_list

            total, outliers, outliers_list = self.summarize_angles()
            mp_data['angles'] = {}
            mp_data['angles']['total'] = total
            mp_data['angles']['outliers'] = outliers
            mp_data['angles']['outliers_list'] = outliers_list

            stats = self.summarize_clashscores()
            mp_data['clashes'] = {}
            mp_data['clashes']['clashscores'] = stats

            total, outliers, outliers_list = self.summarize_clashes()
            mp_data['clashes']['total'] = total
            mp_data['clashes']['outliers'] = outliers
            mp_data['clashes']['outliers_list'] = outliers_list

            stats = self.summarize_rama_rota_stats('rama')
            mp_data['rama'] = {}
            mp_data['rama']['scores'] = stats
            stats = self.summarize_rama_rota('rama')
            mp_data['rama']['outliers_list'] = stats

            stats = self.summarize_rama_rota_stats('rota')
            mp_data['rota'] = {}
            mp_data['rota']['scores'] = stats
            stats = self.summarize_rama_rota('rota')
            mp_data['rota']['outliers_list'] = stats

        return mp_data

    @property
    def models(self) -> list:
        data = self.data
        return list(data.keys())

class MyModelDumper(ihm.dumper._ModelDumper):
    _check = False
    model_id = None

    def __init__(self, model_id: int|None=None):
        self.model_id = model_id

    def dump(self, system, writer):
        seen_types = self.dump_atoms(system, writer)
        self.dump_atom_type(seen_types, system, writer)

    def dump_atoms(self, system, writer, add_ihm=True):

        seen_types = {}
        ordinal = itertools.count(1)
        it = ["group_PDB", "id", "type_symbol", "label_atom_id",
              "label_alt_id", "label_comp_id", "label_seq_id", "auth_seq_id",
              "pdbx_PDB_ins_code", "label_asym_id", "Cartn_x", "Cartn_y",
              "Cartn_z", "occupancy", "label_entity_id", "auth_asym_id",
              "auth_comp_id", "B_iso_or_equiv", "pdbx_PDB_model_num"]
        if add_ihm:
            it.append("ihm_model_id")
        with writer.loop("_atom_site", it) as lp:
            for group, model in system._all_models():
                if self.model_id is not None:
                    if self.model_id != int(model._id):
                        continue
                # rngcheck = _RangeChecker(model, self._check)
                seen_atoms  =  {}

                for atom in model.get_atoms():
                    # rngcheck(atom)
                    seq_id = 1 if atom.seq_id is None else atom.seq_id
                    label_seq_id = atom.seq_id
                    if not atom.asym_unit.entity.is_polymeric():
                        label_seq_id = None
                    comp = atom.asym_unit.sequence[seq_id - 1]
                    seen_types[atom.type_symbol] = None
                    auth_seq_id, ins = \
                        atom.asym_unit._get_auth_seq_id_ins_code(seq_id)

                    # Fix for MolProbity
                    if atom.biso is None:
                        atom.biso = 0.0

                    # Fix for MolProbity
                    if atom.occupancy is None:
                        atom.occupancy = 0.0

                    # Don't write duplicated atoms
                    key = (atom.asym_unit.id, atom.seq_id, atom.atom_id)

                    if key in seen_atoms:
                        logging.warning(f'Skipping duplicated atom {key}')
                        continue
                    else:
                        seen_atoms[key] = None

                    lp.write(id=next(ordinal),
                             type_symbol=atom.type_symbol,
                             group_PDB='HETATM' if atom.het else 'ATOM',
                             label_atom_id=atom.atom_id,
                             label_alt_id=atom.alt_id,
                             label_comp_id=comp.id,
                             label_asym_id=atom.asym_unit._id,
                             label_entity_id=atom.asym_unit.entity._id,
                             label_seq_id=label_seq_id,
                             auth_seq_id=auth_seq_id, auth_comp_id=comp.id,
                             pdbx_PDB_ins_code=ins or ihm.unknown,
                             auth_asym_id=atom.asym_unit.strand_id,
                             Cartn_x=atom.x, Cartn_y=atom.y, Cartn_z=atom.z,
                             B_iso_or_equiv=atom.biso,
                             occupancy=atom.occupancy,
                             pdbx_PDB_model_num=model._id,
                             ihm_model_id=model._id)
        return seen_types


class AtomSiteVariant(ihm.dumper.Variant):
    """Used to select typical PDBx/IHM file output. See :func:`write`."""

    def __init__(self, model_id: int|None=None):
        self.__dumpers = [
            ihm.dumper._EntryDumper()]

        if model_id is not None:
            self.__dumpers.append(MyModelDumper(model_id))
        else:
            self.__dumpers.append(MyModelDumper())

    def get_dumpers(self):
        return self.__dumpers

def _stub():
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="MolProbity validation",
        description="Assess quality of atomic models")

    parser.add_argument('-i', '--input',
                       help='mmCIF file',
                       type=str,
                       required=True)

    parser.add_argument('--cache',
                       help='Cache directory',
                       type=str,
                       default='.')

    parser.add_argument('--nocache',
                       help='Disable cache',
                       action='store_true',
                       default=False)

    args = parser.parse_args()

    I_mp = GetMolprobityInformation(args.input,
                                    cache=args.cache,
                                    nocache=args.nocache)
    d_mp = I_mp.get_mp_data()
