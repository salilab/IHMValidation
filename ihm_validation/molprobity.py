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
from subprocess import run
from mmcif_io import GetInputInformation, MAX_NUM_MODELS
import ihm
import ihm.reader
import collections
import pandas as pd
import csv
import re


class GetMolprobityInformation(GetInputInformation):
    _tempfiles = []

    def __init__(self, mmcif_file, cache):
        super().__init__(mmcif_file)
        self.verify_molprobity_installation()
        self.version = self.get_version()
        self.ID = Path(mmcif_file).stem
        self.nos = min(self.get_number_of_models(), MAX_NUM_MODELS)
        if not Path(cache).is_dir():
            os.makedirs(cache)
            logging.info(f'Created cache directory {cache}')

        self.cache = cache

        self._tempcif = str(Path(self.cache, 'temp.cif'))
        if Path(self._tempcif).is_file():
            os.remove(self._tempcif)
        self.rewrite_mmcif(self._tempcif)
        self._tempfiles.append(self._tempcif)


    def verify_molprobity_installation(self):
        """ Stub for a validation function """
        pass


    def get_version(self, tool: str = 'molprobity.clashscore') -> str:
        """
        Get molprobity version.
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



    def check_for_molprobity(self, filetemp=None) -> bool:
        """ Check the biso and occupancy columns for mmcif files"""
        if filetemp is not None:
            model = ihm.model.Model
            system, = ihm.reader.read(filetemp,
                                      model_class=model)
            models = [
                b for i in system.state_groups for j in i for a in j for b in a]
        else:
            """check if file is in the right format for molprobity analysis """
            models = [
                b for i in self.system.state_groups for j in i for a in j for b in a]
        # and  models[0]._atoms[0].occupancy is not None:
        if models[0]._atoms[0].biso is not None:
            print("File in the appropriate format for molprobity")
            return True
        else:
            print("File is not in the appropriate format for molprobity")
            return False

    def rewrite_mmcif(self, outfn='temp.cif'):
        '''Workaround to generate molprobity-compliant mmCIF'''
        fn = self.mmcif_file
        if Path(outfn).is_file():
            os.remove(outfn)
        encoding = 'utf-8'
        try:
            with open(fn, 'r', encoding=encoding) as f:
                raw = f.readlines()
        except:
            encoding = 'latin-1'
            with open(fn, 'r', encoding=encoding) as f:
                raw = f.readlines()

        loop_lines = []
        loop = False
        out = []
        skip_loop = False
        for line in raw:
            loop_break = False
            if re.match('loop_', line):
                if loop and not skip_loop:
                       out.extend(loop_lines)

                loop = True
                skip_loop = False
                loop_lines = [line]

            elif loop:
                loop_lines.append(line)

                if re.match('_flr', line):
                    skip_loop = True
            else:
                out.append(line)

        if loop and not skip_loop:
                out.extend(loop_lines)

        with open(outfn, 'w', encoding=encoding) as f:
            f.write(''.join(out))



    def check_molprobity_processing(self, output_dict: dict) -> bool:
        """check if molprobity output tables have the same number of lines """
        values = set([len(val) for key, val in output_dict.items()])
        if len(values) == 1:
            return True
        return False

    def run_ramalyze(self, d: dict):
        """run ramalyze to get outliers """
        f_name = str(Path(self.cache, self.ID+'_temp_rama.txt'))
        self._tempfiles.append(f_name)

        with open(f_name, 'w+') as f:
            run(['molprobity.ramalyze', self._tempcif],
                stdout=f,
                cwd=self.cache)

        with open(f_name, 'r') as f:
            line = [_.strip() for _ in f.readlines()]

        d['rama'] = line
        f_name = str(Path(self.cache, self.ID + '_temp_rama.pickle'))

        with open(f_name, 'wb') as f:
            pickle.dump(d['rama'], f)

    def run_molprobity(self, d: dict):
        """run molprobity"""
        f_name = str(Path(
            self.cache, self.ID + '_temp_mp.txt'))
        self._tempfiles.append(f_name)

        with open(f_name, 'w+') as f:
            run(['molprobity.molprobity', self._tempcif,
                 # "disable_uc_volume_vs_n_atoms_check=True",
                 # This is a legacy option and causes extremely
                 # large memory consumption with recent
                 # molprobity versions on PDB-Dev entries
                 "coot=False"],
                stdout=f,
                cwd=self.cache)
            try:
                os.remove(str(Path(self.cache, 'molprobity.out')))
            except OSError:
                logging.error("Couldn't delete molprobity.out")

        with open(f_name, 'r') as f:
            line = [_.strip() for _ in f.readlines()]

        d['molprobity'] = line
        f_name = str(Path(self.cache, self.ID + '_temp_mp.pickle'))

        with open(f_name, 'wb') as f:
            pickle.dump(d['molprobity'], f)

    def run_clashscore(self, d: dict):
        """run clashscore to get information on steric clashes"""
        f_name = str(Path(
            self.cache, self.ID + '_temp_clash.txt'))
        self._tempfiles.append(f_name)

        with open(f_name, 'w+') as f:
            run(['molprobity.clashscore', self._tempcif],
                stdout=f,
                cwd=self.cache)

        with open(f_name, 'r') as f:
            line = [_.strip() for _ in f.readlines()]

        d['clash'] = line

        f_name = str(Path(
            self.cache, self.ID + '_temp_clash.pickle'))

        with open(f_name, 'wb') as f:
            pickle.dump(d['clash'], f)

    def run_rotalyze(self, d: dict):
        """run rotalyZe to get rotameric outliers"""
        f_name = str(self.ID)+'_temp_rota.txt'
        self._tempfiles.append(f_name)

        with open(f_name, 'w+') as f:
            run(['molprobity.rotalyze', self._tempcif],
                stdout=f,
                cwd=self.cache)

        with open(f_name, 'r') as f:
            line = [_.strip() for _ in f.readlines()]

        d['rota'] = line
        f_name = str(Path(self.cache, self.ID+'_temp_rota.pickle'))

        with open(f_name, 'wb') as f:
            pickle.dump(d['rota'], f)

    def write_all_lines(self, file_handle) -> list:
        """print all lines from file to list """
        with open(file_handle.name, 'r') as f:
            line = [_.strip() for _ in f.readlines()]
        return line

    def process_rama(self, line: list) -> dict:
        """ reading and processing molprobity output from rama outliers.
        Outputs information specific to models """
        line_new = line[1:-3]
        count = 1
        models = {_: [] for _ in range(1, self.nos+1)}
        cutoff = len(line_new)/self.nos
        for ind, el in enumerate(line_new):
            if self.nos > 1 and ind < len(line_new)-1:
                if ind < count*cutoff:
                    models[count].append(el)
                else:
                    count = count+1
                    models[count].append(el)
            else:
                models[count].append(el)
        return models

    def process_molprobity(self, line: list) -> (list, list):
        """ process molprobity files to extract relevant information """
        bond_index = angle_index = None
        for ind, el in enumerate(line):
            if 'Bond outliers' in el:
                bond_index = ind
            if 'Angle outliers' in el:
                angle_index = ind
            if 'Molprobity validation' in el:
                end = ind
        if angle_index is None and bond_index is not None:
            bond_outliers = line[bond_index+2:end-1]
            return (bond_outliers, [])
        elif bond_index is None and angle_index is not None:
            angle_outliers = line[angle_index+2:end-1]
            return ([], angle_outliers)
        elif bond_index is None and angle_index is None:
            return ([], [])
        else:
            bond_outliers = line[bond_index+2:angle_index-1]
            angle_outliers = line[angle_index+2:end-1]
            return (bond_outliers, angle_outliers)

    def process_angles(self, line: list) -> (list, int):
        """ process molprobity files to extract relevant information """
        total_angles = 0
        angle_index_beg = None
        angle_index_end = None
        ind_end = len(line)

        def find_end_line(ind_beg, ind_end, match_word):
            for ind in range(ind_beg, ind_end):
                if match_word in line[ind]:
                    return ind-1

        for ind, el in enumerate(line):
            if el.startswith('Bond angle restraints:'):
                total_angles = int(el.split(':', 1)[1])
            if 'Bond angles' in el:
                angle_index_beg = ind+3
                break
        if angle_index_beg:
            angle_index_end = find_end_line(
                angle_index_beg, ind_end, match_word='Min. delta:')
            angle_outliers = line[angle_index_beg:angle_index_end]
            return angle_outliers, total_angles
        else:
            return [], total_angles

    def process_angles_list(self, line: list, chains: list) -> (dict, int):
        """ process molprobity list to dict/table for output """
        angle_outliers, total_angles = self.process_angles(line)
        angle = []

        angledict = {'Number': [], 'Chain': [], 'Residue ID': [],
                     'Residue type': [], 'Angle': [], 'Observed angle (&#176)': [],
                     'Ideal angle (&#176)': [], 'key': [], 'Frequency': []}

        list_for_counter = []

        for ind, outlier in enumerate(angle_outliers):
            sub_line = outlier.split()
            if len(sub_line) == 4 or len(sub_line) == 3:
                angle.append(sub_line[-1])

            elif len(sub_line) == 10:
                angle.append(sub_line[3])
                angledict['Angle'].append('-'.join(angle))
                angledict['Chain'].append(sub_line[0])
                angledict['Residue ID'].append(sub_line[1])
                angledict['Residue type'].append(sub_line[2])
                angledict['Observed angle (&#176)'].append(sub_line[4])
                angledict['Ideal angle (&#176)'].append(sub_line[5])
                angledict['key'].append('-'.join(sub_line[:4]))
                angle = []
                list_for_counter.append('-'.join(sub_line[:4]))

            elif len(sub_line) == 9:
                angle.append(sub_line[2])
                angledict['Angle'].append('-'.join(angle))

                if sub_line[0] in chains:
                    angledict['Chain'].append(sub_line[0])
                    angledict['Residue ID'].append(sub_line[1])
                    angledict['Residue type'].append(sub_line[2])

                elif len(sub_line[0]) > 1 and sub_line[0][:1] in chains:
                    angledict['Chain'].append(sub_line[0][:1])
                    angledict['Residue ID'].append(sub_line[0][1:])
                    angledict['Residue type'].append(sub_line[1])

                elif len(sub_line[0]) > 1 and sub_line[0][:2] in chains:
                    angledict['Chain'].append(sub_line[0][:2])
                    angledict['Residue ID'].append(sub_line[0][2:])
                    angledict['Residue type'].append(sub_line[1])

                angledict['Observed angle (&#176)'].append(sub_line[4])
                angledict['Ideal angle (&#176)'].append(sub_line[5])
                angledict['key'].append('-'.join(sub_line[:4]))
                angle = []
                list_for_counter.append('-'.join(sub_line[:4]))

        freq_dict = collections.Counter(list_for_counter)

        for ind, el in enumerate(angledict['key']):
            angledict['Frequency'].append(freq_dict[el])
            angledict['Number'].append(ind+1)

        del angledict['key']

        if self.check_molprobity_processing(angledict):
            return angledict, total_angles
        else:
            return "Your molprobity processing is incorrect, please check the code", 0

    def add_angles_outliers(self, line: list, angledict: dict, chains: list) -> dict:
        """ add to angle outlier dict/table as molprobity outputs angle outliers in multiple formats """
        list_for_counter = []
        existing_number = len(angledict['Number'])
        angledict['key'] = []
        for ind, outlier in enumerate(line):
            sub_line = outlier.split()

            if sub_line[0] in chains or len(sub_line[0]) == 1:
                angledict['Chain'].append(sub_line[0])
                angledict['Residue ID'].append(sub_line[1])
                angledict['Residue type'].append(sub_line[2])

            else:
                temp = sub_line[0]

                if temp[:1] in chains:
                    val1 = temp[:1]
                    val2 = temp[1:]
                    angledict['Chain'].append(val1)
                    angledict['Residue ID'].append(val2)
                    angledict['Residue type'].append(sub_line[1])

                elif temp[:2] in chains:
                    val1 = temp[:2]
                    val2 = temp[2:]
                    angledict['Chain'].append(val1)
                    angledict['Residue ID'].append(val2)
                    angledict['Residue type'].append(sub_line[1])

            angledict['key'].append('-'.join(sub_line[:4]))

            if sub_line[3] == 'Angle':
                angledict['Angle'].append(sub_line[4][:-1])
                angledict['Observed angle (&#176)'].append(sub_line[6][:-1])
                ideal = round(float(sub_line[6][:-1])+float(sub_line[-1]), 2)
                angledict['Ideal angle (&#176)'].append(ideal)

            elif sub_line[4] == 'Angle':
                angledict['Angle'].append(sub_line[5][:-1])
                angledict['Observed angle (&#176)'].append(sub_line[7][:-1])
                ideal = round(float(sub_line[7][:-1])+float(sub_line[-1]), 2)
                angledict['Ideal angle (&#176)'].append(ideal)

            elif sub_line[5] == 'Angle':
                angledict['Angle'].append(sub_line[6][:-1])
                angledict['Observed angle (&#176)'].append(sub_line[8][:-1])
                ideal = round(float(sub_line[8][:-1])+float(sub_line[-1]), 2)
                angledict['Ideal angle (&#176)'].append(ideal)

            list_for_counter.append('-'.join(sub_line[:4]))

        freq_dict = collections.Counter(list_for_counter)

        for ind, el in enumerate(angledict['key']):
            angledict['Frequency'].append(freq_dict[el])
            angledict['Number'].append(ind+1+existing_number)

        del angledict['key']

        if self.check_molprobity_processing(angledict):
            return angledict
        else:
            return "Your molprobity processing is incorrect, please check the code"

    def angle_summary_table(self, outlier_list: dict) -> dict:
        '''
        converts full detailed list to summary table
        '''
        summary_dict = {'Angle type': [], 'Observed angle (&#176)': [], 'Ideal angle (&#176)': [
        ], 'Number of outliers': []}
        list_for_counter = []

        for ind, val in enumerate(outlier_list[1:]):
            temp = '{0:s}:{1:.2f}:{2:.2f}'.format(
                val[4], float(val[5]), float(val[6]))
            list_for_counter.append(temp)

        freq_dict = collections.Counter(list_for_counter)

        for key, val in freq_dict.items():
            summary_dict['Number of outliers'].append(val)
            summary_dict['Angle type'].append(key.split(':')[0])
            summary_dict['Observed angle (&#176)'].append(key.split(':')[1])
            summary_dict['Ideal angle (&#176)'].append(key.split(':')[2])
        return summary_dict

    def bond_summary_table(self, outlier_list: dict) -> dict:
        '''
        converts full detailed list to summary table
        '''
        summary_dict = {'Bond type': [], 'Observed distance (&#8491)': [
        ], 'Ideal distance (&#8491)': [], 'Number of outliers': []}
        list_for_counter = []

        for ind, val in enumerate(outlier_list[1:]):
            temp = '{0:s}:{1:.2f}:{2:.2f}'.format(
                val[4], float(val[5]), float(val[6]))

            list_for_counter.append(temp)

        freq_dict = collections.Counter(list_for_counter)

        for key, val in freq_dict.items():
            summary_dict['Number of outliers'].append(val)
            summary_dict['Bond type'].append(key.split(':')[0])
            summary_dict['Observed distance (&#8491)'].append(
                key.split(':')[1])
            summary_dict['Ideal distance (&#8491)'].append(key.split(':')[2])
        return summary_dict

    def write_table_csv(self, output_list: list, csvDirName: str, table_filename: str):
        '''
        convert outlier list to csv
        '''
        self.filename = str(Path(csvDirName, table_filename))
        with open(self.filename, 'w') as f:
            write = csv.writer(f)
            for row in output_list:
                write.writerow(row)

    def write_table_html(self, output_list: list, htmlDirName: str, table_filename: str):
        '''
        convert outlier list to html
        '''
        self.filename = str(Path(htmlDirName, table_filename))
        with open(self.filename, 'w') as f:
            f.write('<!DOCTYPE html>\n<html lang="en">\n<body>\n<p>\n')
            for line in output_list:
                f.write(', '.join(line)+'<br>')
            f.write('</p>\n</body>\n</html>')

    def process_bonds(self, line: list) -> (list, int):
        """ process molprobity files to extract relevant information """
        total_bonds = 0
        bond_index_beg = None
        bond_index_end = None
        ind_end = len(line)

        def find_end_line(ind_beg, ind_end, match_word):
            for ind in range(ind_beg, ind_end):
                if match_word in line[ind]:
                    return ind-1

        for ind, el in enumerate(line):
            if el.startswith('Bond restraints:'):
                total_bonds = int(el.split(':', 1)[1])
            if 'Bond outliers' in el:
                bond_index_beg = ind+2
                break
        bond_index_end1 = find_end_line(
            bond_index_beg, ind_end, match_word='Molprobity')
        bond_index_end2 = find_end_line(
            bond_index_beg, ind_end, match_word='Angle')
        bond_index_end = min(bond_index_end1, bond_index_end2)
        bond_outliers = line[bond_index_beg:bond_index_end]
        return bond_outliers, total_bonds

    def process_bonds_list(self, line: list, chains: list) -> (dict, int):
        """ process molprobity files to extract relevant information """

        bond_outliers, total_bonds = self.process_bonds(line)
        bonddict = {'Number': [], 'Chain': [], 'Residue ID': [],
                    'Residue type': [], 'Bond': [], 'Observed distance (&#8491)': [],
                    'Ideal distance (&#8491)': [], 'key': [], 'Frequency': []}

        list_for_counter = []

        for ind, outlier in enumerate(bond_outliers):
            sub_line = outlier.split()

            if sub_line[0] in chains or len(sub_line[0]) == 1:
                bonddict['Chain'].append(sub_line[0])
                bonddict['Residue ID'].append(sub_line[1])
                bonddict['Residue type'].append(sub_line[2])

            else:
                temp = sub_line[0]

                if temp[:1] in chains:
                    val1 = temp[:1]
                    val2 = temp[1:]
                    bonddict['Chain'].append(val1)
                    bonddict['Residue ID'].append(val2)
                    bonddict['Residue type'].append(sub_line[1])

                elif temp[:2] in chains:
                    val1 = temp[:2]
                    val2 = temp[2:]
                    bonddict['Chain'].append(val1)
                    bonddict['Residue ID'].append(val2)
                    bonddict['Residue type'].append(sub_line[1])

            bonddict['key'].append('-'.join(sub_line[:4]))
            bonddict['Number'].append(ind+1)

            if sub_line[3] == 'Bond':
                bonddict['Bond'].append(sub_line[4][:-1])
                bonddict['Observed distance (&#8491)'].append(sub_line[6][:-1])
                ideal = round(float(sub_line[6][:-1])+float(sub_line[-1]), 2)
                bonddict['Ideal distance (&#8491)'].append(ideal)

            elif sub_line[4] == 'Bond':
                bonddict['Bond'].append(sub_line[5][:-1])
                bonddict['Observed distance (&#8491)'].append(sub_line[7][:-1])
                ideal = round(float(sub_line[7][:-1])+float(sub_line[-1]), 2)
                bonddict['Ideal distance (&#8491)'].append(ideal)

            elif sub_line[5] == 'Bond':
                bonddict['Bond'].append(sub_line[6][:-1])
                bonddict['Observed distance (&#8491)'].append(sub_line[8][:-1])
                ideal = round(float(sub_line[8][:-1])+float(sub_line[-1]), 2)
                bonddict['Ideal distance (&#8491)'].append(ideal)

            list_for_counter.append('-'.join(sub_line[:4]))

        freq_dict = collections.Counter(list_for_counter)

        for el in bonddict['key']:
            bonddict['Frequency'].append(freq_dict[el])

        del bonddict['key']

        if self.check_molprobity_processing(bonddict):
            return bonddict, total_bonds
        else:
            return "Your molprobity processing is incorrect, please check the code", 0

    @staticmethod
    def get_model_id_str(line: str) -> str:
        """ extract MODEL X substring """
        m = re.search('MODEL\s*(?P<model_id>\d+)', line, re.IGNORECASE)

        mid = None

        if m:
            g = m.groupdict()
            mid = g['model_id']

        return mid


    def process_clash(self, line: list) -> dict:
        """ process clash files to extract relevant information """
        count = [i for i, j in enumerate(line) if 'Bad Clashes' in j]
        if self.nos > 1:
            clashes = {f'Model {self.get_model_id_str(j)}':[]
                    for i, j in enumerate(line) if 'Bad Clashes' in j}
        else:
            clashes = {'Model 1': []}
        count.append(self.find_clashscore_records(line))

        for ind in range(0, len(clashes.keys())):
            output_line = [j for k, j in enumerate(line) if k > int(
                count[ind]) and k < int(count[ind+1])]
            clashes[list(clashes.keys())[ind]].append(output_line)

        return clashes

    @staticmethod
    def find_clashscore_records(line: list) -> int:
        """ look for the beginning of clashscore records in clash data.
        searches for strings like 'MODEL 1 clashscore = 42.58'
        """
        found = False
        q = re.compile('clashscore')
        for i, line_ in enumerate(line):
            if q.search(line_):
                found = True
                break

        if not found:
            logging.warning('Could not find clashscore records')

        return i

    def process_rota(self, line: list) -> dict:
        """ process rota files to extract relevant information """
        line_new = line[1:-1]
        count = 1
        models = {_: [] for _ in range(1, self.nos+1)}
        cutoff = len(line_new)/self.nos
        for ind, el in enumerate(line_new):
            if self.nos > 1 and ind < len(line_new)-1:
                if ind < count*cutoff:
                    models[count].append(el)
                else:
                    count = count+1
                    models[count].append(el)
            else:
                models[count].append(el)
        return models

    def rama_summary_table(self, models: dict) -> dict:
        """ write out summary table from rama, clash and other tables"""
        f_rama = open(str(Path(self.cache, self.ID+'_rama_summary.txt')), 'w+')
        dict1 = {'Model ID': [], 'Analyzed': [],
                 'Favored': [], 'Allowed': [], 'Outliers': []}
        for ind, el in models.items():
            dict1['Model ID'].append(ind)
            F = []
            A = []
            U = []
            for line in el:
                if line.strip().split()[-2] == 'or':
                    subline = ':'.join(line.strip().split()[-3:])
                else:
                    subline = line.strip().split()[-1]
                if subline.split(':')[4] == 'Favored':
                    F.append(subline.split(':')[4])
                elif subline.split(':')[4] == 'Allowed':
                    A.append(subline.split(':')[4])
                else:
                    U.append(subline.split(':')[4])
            dict1['Analyzed'].append(len(F)+len(A)+len(U))
            dict1['Favored'].append(len(F))
            dict1['Allowed'].append(len(A))
            dict1['Outliers'].append(len(U))
        print(dict1['Model ID'], file=f_rama)
        print(dict1['Analyzed'], file=f_rama)
        print(dict1['Favored'], file=f_rama)
        print(dict1['Allowed'], file=f_rama)
        print(dict1['Outliers'], file=f_rama)
        return dict1

    def clash_summary_table(self, line: list) -> (dict, int):
        """ format clash data to print to file """
        def get_clash_score(line: str) -> str:
            """ parse clash line """
            m = re.search(
                'clashscore\s+=\s+(?P<clashscore>\d+\.\d+)',
                line,
                re.IGNORECASE
            )

            g = m.groupdict()

            return g['clashscore']

        with open(
            str(Path(self.cache,  self.ID+'_clash_summary.txt')), 'w+') as f_clash:

            clashes = self.process_clash(line)
            if self.nos > 1:
                # Find the beginning of clashscore records
                cs_start = self.find_clashscore_records(line)
                # Extract only X models
                clashscore_list = line[cs_start:cs_start + self.nos]
            else:
                clashscore_list = ['Model 1 ' + (line[len(line)-self.nos:])[0]]
            dict1 = {'Model ID': [], 'Clash score': [], 'Number of clashes': []}

            for clashval in clashscore_list:
                mid = self.get_model_id_str(clashval)
                clashscore = get_clash_score(clashval)
                dict1['Model ID'].append(mid)
                dict1['Clash score'].append(clashscore)

            for model_id in dict1['Model ID']:
                dict1['Number of clashes'].append(len(clashes[f'Model {model_id}'][0]))
            clash_total = (sum(dict1['Number of clashes']))
            print(dict1['Model ID'], file=f_clash)
            print(dict1['Clash score'], file=f_clash)
            print(dict1['Number of clashes'], file=f_clash)

        return dict1, clash_total

    def orderclashdict(self, modeldict: dict) -> dict:
        """ DEPRECATED molprobity returns output in lexicographic order
        this is to change it to number order
         """
        df = pd.DataFrame(modeldict)
        df['ID'] = df['Model ID'].apply(lambda x: int(x.split()[1]))
        df = df.sort_values(by='ID')
        df = df.drop(['ID'], axis=1)
        df_dict = df.to_dict()
        ordered_dict = {key: list(val.values())
                        for key, val in df_dict.items()}
        return ordered_dict

    def rota_summary_table(self, models: dict) -> dict:
        """ format rota data to print to file """
        dict1 = {'Model ID': [], 'Analyzed': [],
                 'Favored': [], 'Allowed': [], 'Outliers': []}
        f_rota = open(str(Path(self.cache, self.ID+'_rota_summary.txt')), 'w+')
        for ind, el in models.items():
            dict1['Model ID'].append(ind)
            F = []
            A = []
            U = []
            for line in el:
                if line.strip().split()[-1].split(':')[-2] == 'Favored':
                    F.append(line.strip().split()[-1].split(':')[-2])
                elif line.strip().split()[-1].split(':')[-2] == 'Allowed':
                    A.append(line.strip().split()[-1].split(':')[-2])
                else:
                    U.append(line.strip().split()[-1].split(':')[-2])
            dict1['Analyzed'].append(len(F)+len(A)+len(U))
            dict1['Favored'].append(len(F))
            dict1['Allowed'].append(len(A))
            dict1['Outliers'].append(len(U))
        print(dict1['Model ID'], file=f_rota)
        print(dict1['Analyzed'], file=f_rota)
        print(dict1['Favored'], file=f_rota)
        print(dict1['Allowed'], file=f_rota)
        print(dict1['Outliers'], file=f_rota)
        return dict1

    def rama_detailed_table(self, models: dict, chains: list) -> dict:
        """ format rama information to print to file"""
        dict1 = {'Model ID': [], 'Chain': [],
                 'Residue ID': [], 'Residue type': []}
        f_rama_D = open(str(Path(self.cache, self.ID+'_rama_detail.txt')), 'w+')
        for ind, el in models.items():
            for line in el:
                if line.strip().split()[-2] == 'or':
                    subline = ':'.join(line.strip().split()[-3:])
                else:
                    subline = line.strip().split()[-1]
                if subline.split(':')[4] == 'OUTLIER':
                    dict1['Model ID'].append(ind)
                    dict1['Residue type'].append(subline.split(':')[0])
                    if len(line.strip().split()[0]) > 2:
                        temp = line.strip().split()[0]
                        if temp[:1] in chains:
                            val1 = temp[:1]
                            val2 = temp[1:]
                            dict1['Chain'].append(val1)
                            dict1['Residue ID'].append(val2)
                        elif temp[:2] in chains:
                            val1 = temp[:2]
                            val2 = temp[2:]
                            dict1['Chain'].append(val1)
                            dict1['Residue ID'].append(val2)
                        elif temp[:3] in chains:
                            val1 = temp[:3]
                            val2 = temp[3:]
                            dict1['Chain'].append(val1)
                            dict1['Residue ID'].append(val2)
                    else:
                        val1 = line.strip().split()[0]
                        val2 = line.strip().split()[1]
                        dict1['Chain'].append(val1)
                        dict1['Residue ID'].append(val2)

        print(dict1['Model ID'], file=f_rama_D)
        print(dict1['Chain'], file=f_rama_D)
        print(dict1['Residue ID'], file=f_rama_D)
        print(dict1['Residue type'], file=f_rama_D)
        return dict1

    def clash_detailed_table(self, line: list) -> dict:
        """process molprobity clash information and format to table"""
        f_clash_D = open(str(Path(self.cache, self.ID+'_clash_detailed.txt')), 'w+')
        dict1 = {'Model ID': [], 'Atom-1': [],
                 'Atom-2': [], 'Clash overlap (&#8491)': []}
        clashes = self.process_clash(line)
        for ind, el in clashes.items():
            for line in el[0]:
                subline = [_ for _ in line.split(' ') if _ not in '']
                dict1['Model ID'].append(self.get_model_id_str(ind))
                if len(subline) < 9 and len(subline[0]) > 2:
                    dict1['Atom-1'].append(':'.join(subline[0:3]))
                else:
                    dict1['Atom-1'].append(':'.join(subline[0:4]))
                if len(subline) < 9 and len(subline[3]) > 4:
                    dict1['Atom-2'].append(':'.join(subline[3:-1]))
                else:
                    dict1['Atom-2'].append(':'.join(subline[4:-1]))
                dict1['Clash overlap (&#8491)'].append(
                    subline[-1].replace(':', ''))

        print(dict1['Model ID'], file=f_clash_D)
        print(dict1['Atom-1'], file=f_clash_D)
        print(dict1['Atom-2'], file=f_clash_D)
        print(dict1['Clash overlap (&#8491)'], file=f_clash_D)
        return dict1

    def rota_detailed_table(self, models: dict, chains: list) -> dict:
        """process molprobity rotamers information and format to table"""
        f_rota_D = open(str(Path(self.cache,
                                     self.ID+'_rota_detailed.txt')), 'w+')
        dict1 = {'Model ID': [], 'Chain': [],
                 'Residue ID': [], 'Residue type': []}
        for ind, el in models.items():
            for line in el:
                if line.strip().split()[-1].split(':')[-2] == 'OUTLIER':
                    dict1['Model ID'].append(ind)
                    dict1['Residue type'].append(
                        line.strip().split()[-1].split(':')[0])
                    if len(line.strip().split()[0]) > 2:
                        temp = line.strip().split()[0]
                        if temp[:1] in chains:
                            val1 = temp[:1]
                            val2 = temp[1:]
                            dict1['Chain'].append(val1)
                            dict1['Residue ID'].append(val2)

                        elif temp[:2] in chains:
                            val1 = temp[:2]
                            val2 = temp[2:]
                            dict1['Chain'].append(val1)
                            dict1['Residue ID'].append(val2)

                        elif temp[:3] in chains:
                            val1 = temp[:3]
                            val2 = temp[3:]
                            dict1['Chain'].append(val1)
                            dict1['Residue ID'].append(val2)

                    else:
                        val1 = line.strip().split()[0]
                        val2 = line.strip().split()[1]
                        dict1['Chain'].append(val1)
                        dict1['Residue ID'].append(val2)

        print(dict1['Model ID'], file=f_rota_D)
        print(dict1['Chain'], file=f_rota_D)
        print(dict1['Residue ID'], file=f_rota_D)
        print(dict1['Residue type'], file=f_rota_D)
        return dict1

    def get_data_for_quality_at_glance(self, clash: list, rota: list, rama: list) -> dict:
        """format mean information of models for quality at glance plots, read from temp_mp file"""
        molprobity = {'Names': [], 'Models': [], 'Clashscore': [],
                      'Ramachandran outliers': [], 'Sidechain outliers': []}
        for ind, model in enumerate(clash):
            if ind > 0:
                molprobity['Names'].append(model[0])
                molprobity['Models'].append(ind)
                if len(model[1]) > 0:
                    molprobity['Clashscore'].append(round(float(model[1]), 2))
                else:
                    molprobity['Clashscore'].append(0.0)
                molprobity['Ramachandran outliers'].append(int(rama[ind][-1]))
                molprobity['Sidechain outliers'].append(int(rota[ind][-1]))
        return molprobity

    def cleanup(self):
        for f_name in self._tempfiles:
            try:
                os.remove(f_name)
            except OSError:
                logging.error(f"Coldn't delete temp file: {f_name}")
