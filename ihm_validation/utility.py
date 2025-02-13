###################################
# Script:
# 1) Contains trivial and useful
# functions
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

import os
import glob
from pathlib import Path
from collections import Counter, defaultdict
from multiprocessing import Process
import numpy as np
import logging
import ihm, ihm.reader, ihm.model
import itertools
import time
import signal
import re
import requests

NA = 'Not available'

def dict_to_JSlist(d: dict) -> list:
    '''
    convert dictionary to list of lists
    '''
    output_list = []

    if bool(d) and len(list(d.keys())) > 0:
        # add headers for table, which are the keys of the dict
        header = list(d.keys())
        # Get number of columns
        N = len(header)
        # Get number of rows. +1 for header.
        M = len(d[header[0]]) + 1
        output_list = np.empty((M, N), dtype=object)
        output_list[0, :] = header
        # iterate over dict keys - columns
        for j, v in enumerate(d.values()):
            # iterate over values of every key - fill rows
            for i, el in enumerate(v, start=1):
                # Check if int or float
                if isinstance(el, int) or isinstance(el, float):
                    el_ = el

                # If string, try casting as int or float
                elif isinstance(el, str):
                    try:
                        el_ = int(el)
                    except (TypeError, ValueError):
                        el_ = str(el)

                    if isinstance(el_, str):
                        try:
                            el_ = float(el)
                        except (TypeError, ValueError):
                            el_ = str(el)

                # Otherwise cast as str
                else:
                    el_ = str(el)

                if el_ == '?':
                    el_ = '_'
                try:
                    output_list[i, j] = el_
                except IndexError:
                    logging.error(
                        'Dict has excessive elements. Ignoring them.')

        output_list = output_list.tolist()

    return output_list


def format_RB_text(tex: list) -> str:
    '''
    convert RB information to text for supp table
    '''
    val = []
    for el in tex:
        for subel in el:
            if subel == el[-1] and el == tex[-1]:
                val.append(str(subel))
            elif subel == el[-1] and el != tex[-1]:
                val.append(str(subel)+', ')
            else:
                val.append(str(subel)+':')

    if len(val) == 0 or len(tex) == 0 or val == '':
        val = ['-']
    return ''.join(val)


def format_flex_text(tex: list) -> str:
    '''
    convert flex information to text for supp table
    '''
    val = []
    for el in tex:
        for subel in el:
            if subel == el[-1] and el == tex[-1]:
                val.append(str(subel))
            else:
                val.append(str(subel)+', ')

    if len(val) == 0 or len(tex) == 0 or val == '':
        val = ['-']

    return ''.join(val)


def format_tuple(tex: list) -> str:
    return str(tex[0])+'-'+str(tex[1])


def dict_to_JSlist_rows(dict1: dict, dict2: dict) -> list:
    '''
    format rigid and flexible segments
    '''
    output_list = []
    output_list.append(['Chain ID', 'Rigid segments', 'Flexible segments'])
    for ind, el in dict1.items():
        output_list.append(
            [ind, format_RB_text(el), format_flex_text(dict2[ind])])
    return output_list


def islistempty(inlist: list) -> bool:
    '''
    minor func
    '''
    if isinstance(inlist, list):
        return all(map(islistempty, inlist))
    return False


def cat_list_string(listn: list) -> str:
    '''
    minor func
    '''
    result = ' '
    for ind in range(len(listn)):
        if ind == 0:
            result += str(listn[ind])
        else:
            result += ','
            result += str(listn[ind])
    return result


def get_key_from_val(dict1: dict, val1: str) -> list:
    '''
    minor func
    '''
    return dict1.keys()[dict1.values().index(val1)]


def get_val_from_key(dict1: dict, key1: str) -> list:
    '''
    minor func
    '''
    return dict1[key1]


def get_name(name) -> str:
    '''
    minor func
    '''
    return str(name)


def get_copy(name):
    '''
    minor func
    '''
    if str(name) == '?':
        copy = 'None listed'
    elif '.' in name:
        copy = (name.split('.')[1]).split('@')[0]
    else:
        copy = 1
    return copy


def get_unique_datasets(name: dict) -> list:
    '''
    get all datatypes that are yet to be validated
    the ones that can't or the ones that have already been validated
    are in the sub_data set
    '''
    all_data = set(name['Dataset type'])
    sub_data = {'Integrative model', 'Other', 'Comparative model',
                'Experimental model', 'De Novo model', 'SAS data', 'Crosslinking-MS data'}
    fin_data = list(all_data.difference(sub_data))
    output = list()
    for i in fin_data:
        min_list = [j for j in i.split() if j not in ['data']]
        output.append(' '.join(min_list))
    return output


def get_all_files(path_dir):
    '''
    minor func
    '''
    return glob.glob(path_dir)


def runInParallel(*fns):
    '''
    minor func
    '''
    proc = []
    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()


def runInParallel_noargs(*fns):
    '''
    minor func
    '''
    proc = []
    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()


def get_output_file_html(prefix: str) -> str:
    '''
    minor func
    '''
    return f'ValidationReport_{prefix}.html'


def get_supp_file_html(prefix: str) -> str:
    '''
    minor func
    '''
    return f'Supplementary_{prefix}.html'


def get_output_file_temp_html(prefix: str) -> str:
    '''
    minor func
    '''
    return 'temp.html'


def get_output_file_pdf(prefix: str) -> str:
    '''
    minor func
    '''
    return f'{prefix}_full_validation.pdf'


def get_output_file_json(prefix: str) -> str:
    '''
    minor func
    '''
    return f'ValidationReport_{prefix}.json'


def get_supp_file_pdf(prefix: str) -> str:
    '''
    minor func
    '''
    return f'{prefix}_summary_validation.pdf'


def get_subunits(sub_dict: dict) -> list:
    '''
    format chains for supplementary/summary table
    '''
    model_number = len(sub_dict['Model ID'])
    sublist = ['%s: Chain %s (%s residues)' % (sub_dict['Subunit name'][i], sub_dict['Chain ID']
                                               [i], str(sub_dict['Total residues'][i])) for i in range(model_number)]
    return list(set(sublist))


def get_datasets(data_dict: dict) -> list:
    '''
    format datasets for supplementary/summary table
    '''

    dataset_number = len(data_dict['ID'])
    datalist = ['%s, %s' % (data_dict['Dataset type'][i], data_dict['Details'][i])
                for i in range(dataset_number)]
    return datalist


def get_software(data_dict: dict) -> list:
    '''
    format software for supplementary/summary table
    '''

    if len(data_dict) > 0:
        dataset_number = len(data_dict['ID'])
        datalist = ['%s (version %s)' % (data_dict['Software name'][i],
                                         data_dict['Software version'][i]) for i in range(dataset_number)]
        return datalist
    else:
        return ['Software details not provided']


def get_RB(data_list: list) -> list:
    '''
    format RB for supplementary/summary table
    '''

    data_num = len(data_list)
    # datalist = ['%s: %s ' % (data_list[i][0], data_list[i][1],)
    #           for i in range(1, data_num)]
    datalist = []
    for i in range(1, data_num):
        if len(data_list[i][1]) < 1:
            data_list[i][1] = 'None'
        datalist.append('%s: %s ' % (data_list[i][0], data_list[i][1]))
    return datalist


def get_flex(data_list: list) -> list:
    '''
    format flexible regions for supplementary/summary table
    '''

    data_num = len(data_list)
    datalist = ['%s: %s ' % (data_list[i][0], data_list[i][2],)
                for i in range(1, data_num)]
    return datalist


def get_method_name(sample_dict: dict) -> str:
    '''
    format method name for supplementary/summary table
    '''

    datastr = '%s ' % (sample_dict['Method name'][0])
    return datastr.replace('monte carlo', 'Monte Carlo')


def get_method_type(sample_dict: dict) -> str:
    '''
    format method type  for supplementary/summary table
    '''

    datastr = '%s ' % (sample_dict['Method type'][0])
    return datastr.replace('monte carlo', 'Monte Carlo')


def get_restraints_info(restraints: dict) -> list:
    '''
    format restraints info for supplementary/summary table
    '''

    restraints_num = len(restraints['Restraint type'])
    datalist = []
    try:
        dataset = [(restraints['Restraint info'][i], restraints['Restraint type'][i])
                   for i in range(restraints_num)]
    except (ValueError, TypeError, IndexError):
        new_restraints = {key: list(set(val))
                          for key, val in restraints.items()}
        restraints_num = min(len(new_restraints['Restraint info']), len(
            new_restraints['Restraint type']))
        dataset = [(new_restraints['Restraint info'][i], new_restraints['Restraint type'][i])
                   for i in range(restraints_num)]
    for i, j in Counter(dataset).items():
        datalist.append(['%s unique %s: %s' % (j, i[1], i[0])])
    return datalist


def format_list_text(sublist: list) -> str:
    '''
    minor func
    '''
    val = ''
    for el in sublist:
        if el == sublist[-1]:
            val += str(el)+'. '
        else:
            val += str(el)+', '
    if val == '':
        val = '-'
    return val


def all_same(items: list):
    '''
    minor func
    '''
    return all(x == items[0] for x in items)


def mp_readable_format(mp: dict) -> list:
    '''
    Format MolProbity results for supplementary/summary table
    '''
    fin_string = []
    for ind, el in enumerate(mp['Models']):
        fin_string.append('Model-'+str(el)+': '+'Clashscore = ' +
                          str(mp['Clashscore'][ind]) + ', ' + 'Number of Ramachandran outliers = ' +
                          str(mp['Ramachandran outliers'][ind]) + ', '+'Number of sidechain outliers = ' +
                          str(mp['Sidechain outliers'][ind]))
    return fin_string


def get_rg_data(rg_dict: dict) -> list:
    '''
    format rg data for supplementary/summary table
    '''

    fin_rg = []
    for key, val in rg_dict.items():
        fin_rg.append(key+': Rg from Gunier is ' +
                      str(val[0])+'nm and Rg from p(r) is ' + str(val[1])+'nm')
    return fin_rg


def get_rg_data_fits(rg_dict: dict) -> list:
    '''
    format sas model fits for supplementary/summary table
    '''

    fin_rg = []
    for key, val in rg_dict.items():
        for ind, el in enumerate(val):
            count = ind+1
            fin_rg.append(key+': Fit ' + str(count) +
                          ' with &#x3A7;&#xb2; value ' + str(el))
    return fin_rg


def get_cx_data_fits(cx_dict: dict) -> list:
    '''
    format crosslinking-MS data for supplementary/summary table
    '''

    fin_cx = []
    count = 0
    for key, val in cx_dict.items():
        count += 1
        fin_cx.append('Crosslinking-MS Fit of medioid: model # ' + str(count) +
                      ', percentage satisfied ' + str(round(val, 2))+'%')
    return fin_cx


def clean_all(report=None):
    '''
    delete all generated files
    '''

    # dirname_ed = os.getcwd()
    # os.listdir('.')
    # for item in os.listdir('.'):
    #     if item.endswith('.txt'):
    #        os.remove(item)
    #    if item.endswith('.csv'):
    #        os.remove(item)
    #    if item.endswith('.json'):
    #        os.remove(item)
    #    if item.endswith('.sascif'):
    #        os.remove(item)

    if report:
        report.clean()


# Adapted from
# https://stackoverflow.com/questions/52859751/
# most-efficient-way-to-find-order-of-magnitude-of-float-in-python
def order_of_magnitude(value: float) -> float:
    '''
    calculate the order of magnitude for a given number

    >>> order_of_magnitude(135)
    2.0

    '''
    if value <= 0:
        raise(f'Wrong value: {value}. '
              'This function works only for positive values')
    return np.floor(np.log10(value))


def calc_optimal_range(counts: list) -> tuple:
    '''
    heuristics to find optimal range for plots

    >>> calc_optimal_range((10, 1567))
    (9.0, 1568.5669999999998)

    '''

    # Find min/max values
    upper = max(counts)
    lower = min(counts)

    # In peculiar cases add an arbitary offset to the range
    if upper == 0:
        upper = 10
        lower = 0

    # Find the data's order of magnitude and make offset.
    # Typically it would be .001%
    if upper > 0:
        oom = order_of_magnitude(upper)
        upper = upper * (1 + 10 ** (-oom))

    # if lower > 0:
    #     oom = order_of_magnitude(lower)
    #     # Do not allow the range to go below zero
    #     lower = max(0, lower * (1 - 10 ** (-oom)))

    lower = 0

    assert lower >= 0 and upper > 0

    return(lower, upper)

def compress_cx_stats(cx_stats: dict) -> list:
    '''Extract per-model satisfactions stats as a flat list'''
    out_stats = []
    for sg, sgv in cx_stats.items():
        for st, stv in sgv.items():
            for mg, mgv in stv.items():
                out_stats.append(mgv['cx_stats']['All']['Satisfied'])

    return out_stats

def get_python_ihm_version() -> str:
    """returns Python-IHM version"""
    import ihm
    return ihm.__version__

def get_hierarchy_from_atoms(atoms) -> dict:
    """Construct polymer hierarchy from a list of atoms"""
    def infinite_defaultdict(): return defaultdict(infinite_defaultdict)
    root = infinite_defaultdict()

    for a in atoms:
        root[a.asym_unit.id][a.seq_id][a.atom_id] = a

    return root

def get_hierarchy_from_model(model) -> dict:
    """Construct polymer hierarchy from atoms and beads in the model"""
    def infinite_defaultdict(): return defaultdict(infinite_defaultdict)
    root = infinite_defaultdict()

    for a in model.get_atoms():
        root[a.asym_unit.id][a.seq_id][a.atom_id] = a

    for r in model.representation:
        if r.granularity == 'by-residue':
            for i in range(r.asym_unit.seq_id_range[0],
                           r.asym_unit.seq_id_range[1] + 1):
                root[r.asym_unit.asym.id][i]['CA'] = None
                root[r.asym_unit.asym.id][i]['coarse-grained'] = None

        elif r.granularity == 'by-feature':
            for i in range(r.asym_unit.seq_id_range[0],
                           r.asym_unit.seq_id_range[1] + 1):
                root[r.asym_unit.asym.id][i]['coarse-grained'] = None

    for s in model.get_spheres():
        # Consider only by-residue spheres
        bs = get_bead_size(s)
        if bs == 1:

            seq_id = s.seq_id_range[0]

            if root[s.asym_unit.id][seq_id]['CA'] is None:
                root[s.asym_unit.id][seq_id]['CA'] = s
            if root[s.asym_unit.id][seq_id]['coarse-grained'] is None:
                root[s.asym_unit.id][seq_id]['coarse-grained'] = s

        else:

            for seq_id in range(s.seq_id_range[0], s.seq_id_range[1] + 1):

                if root[s.asym_unit.id][seq_id]['coarse-grained'] is None:
                    root[s.asym_unit.id][seq_id]['coarse-grained'] = s
                else:
                    s_ = root[s.asym_unit.id][seq_id]['coarse-grained']
                    # Select best possible resolution
                    if get_bead_size(s) < get_bead_size(s_):
                        root[s.asym_unit.id][seq_id]['coarse-grained'] = s


    return root

def get_bead_size(sphere: ihm.model.Sphere) -> int:
    """Number of residues per bead"""
    return sphere.seq_id_range[1] - sphere.seq_id_range[0] + 1

def pretty_print_representations(reprs: dict) -> list:
    """Pretty print information about representation scales"""
    pretty_reprs = []
    for reprs_ in reprs:
        out = ''
        if (reprs_['atomic'] and reprs_['coarse-grained']) or \
        (reprs_['coarse-grained'] and len(reprs_['coarse-grain_levels']) > 1):
            out += 'Multiscale: '

        if reprs_['atomic']:
            out += 'Atomic'

        if reprs_['coarse-grained']:
            if out != '':
                if out[-1] != ' ':
                    out += '; '

            out += 'Coarse-grained: '
            levels = reprs_['coarse-grain_levels']
            if len(levels) == 1:
                out += f'{levels[0]:d}'
            else:
                min_level = min(levels)
                max_level = max(levels)
                out += f'{min_level:d} - {max_level:d}'
            out += ' residue(s) per bead'

        pretty_reprs.append(out)

    return pretty_reprs

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]

def check_for_dataset_type(dataset_list: list=None, dataset_type=None) -> bool:
    """check if the specific dataset type is present in the dataset list"""
    flag = False

    for dataset in dataset_list:
        if isinstance(dataset, dataset_type):
            flag = True

    return flag


def summarize_entities(rep_info: dict) -> list:
    sum_entities = []

    for rep_ in rep_info:
        for k, v in rep_['Chains'].items():
            if isinstance(v['Total residues'], int):
                tr = f'{v["Total residues"]} residues'
            else:
                tr = None
            data_ = (v['Molecule name'], ', '.join(v['Chains']), tr)
            sum_entities.append(data_)

    sum_entities = sorted(set(sum_entities), key=lambda x: x[1])

    output = []

    for e in sum_entities:
        if e[2] is None:
            l = f"{e[0]}: chain(s) {e[1]}"
        else:
            l = f"{e[0]}: chain(s) {e[1]} ({e[2]})"

        output.append(l)

    return output

def summarize_segments(rep_info: dict) -> list:
    output = []

    for rep_ in rep_info:
        rigid, flexible = 0, 0
        for k, v in rep_['Chains'].items():
            rigid_ = len(v['Rigid segments']) * len(v['Chains'])
            flexible_ = len(v['Flexible segments']) * len(v['Chains'])

            rigid += rigid_
            flexible += flexible_

        output.append((f'{rigid}, {flexible}'))

    return output

def parse_ihm_cif(fname, encoding='utf8') -> tuple:
    try:
        with open(fname, encoding=encoding) as fh:
            system, = ihm.reader.read(fh)
    except UnicodeDecodeError:
        encoding = 'ascii'
        with open(fname, encoding=encoding, errors='ignore') as fh:
            system, = ihm.reader.read(fh)

    return(system, encoding)

def is_atomic(data: ihm.System|ihm.model.Model):
    flag = False

    if isinstance(data, ihm.System):
        for group, model in data._all_models():
            flag = is_atomic(model)
            if flag:
                break

    elif isinstance(data, ihm.model.Model):
        if len(data._atoms) > 0:
            flag = True

    return flag

def is_cg(data: ihm.System|ihm.model.Model):
    flag = False

    if isinstance(data, ihm.System):
        for group, model in data._all_models():
            flag = is_cg(model)
            if flag:
                break

    elif isinstance(data, ihm.model.Model):
        if len(data._spheres) > 0:
            flag = True

    return flag

# https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish
class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def format_range(data, format=".2f"):
    if len(data) > 1:
        min_ = np.nanmin(data)
        max_ = np.nanmax(data)
        r_ = f'{min_:{format}}-{max_:{format}}'
    else:
        min_ = np.nanmin(data)
        r_ = f'{min_:{format}}'

    return r_

def get_alphafolddb_link(acc: str) -> str|None:
    """Format link for AlphaFold DB"""
    url = None
    regexp = '^AF-(?P<uniprot>[0-9A-Za-z]+)-F1$'
    m = re.match(regexp, acc)
    if m:
        uid = m.groupdict()['uniprot']
        r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}")
        if r.status_code == 200:
            url = f"https://alphafold.ebi.ac.uk/entry/{uid}"

    return url
