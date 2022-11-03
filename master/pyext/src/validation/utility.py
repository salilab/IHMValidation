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
from collections import Counter
from multiprocessing import Process
import numpy as np
import logging


def dict_to_JSlist_v0(d: dict) -> list:
    '''
    TO BE REMOVED
    convert dictionary to list of lists
    output_list = []
    if bool(d):
        if len(list(d.keys()))>0:
            output_list.append(list(d.keys()))
            target=list(d.values())
            for ind in range(len(target[0])):
                sublist=[]
                for el in target:
                    sublist.append(str(el[ind]))
                output_list.append(sublist)
    return output_list
    '''
    output_list = []
    if bool(d) and len(list(d.keys())) > 0:
        # add headers for table, which are the keys of the dict
        output_list.append(list(d.keys()))
        # add each row of the table as a list
        target = list(d.values())
        for ind in range(len(target[0])):
            sublist = []
            for el in target:
                el = ['_' if str(i) == '?' else str(i) for i in el]
                sublist.append(str(el[ind]))
            output_list.append(sublist)
    return output_list


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
                val.append(str(subel)+'. ')
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
                val.append(str(subel)+'. ')
            else:
                val.append(str(subel)+', ')

    if len(val) == 0 or len(tex) == 0 or val == '':
        val = ['-']

    return ''.join(val)


def format_tuple(tex: list) -> str:
    return str(tex[0])+'-'+str(tex[1])


def dict_to_JSlist_rows(dict1: dict, dict2: dict) -> list:
    '''
    format rigid bodies and flexible elements
    '''
    output_list = []
    output_list.append(['Chain ID', 'Rigid bodies', 'Non-rigid segments'])
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
                'Experimental model', 'De Novo model', 'SAS data'}
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


def get_output_file_html(mmcif_file: str) -> str:
    '''
    minor func
    '''
    return 'ValidationReport_'+mmcif_file+'.html'


def get_supp_file_html(mmcif_file: str) -> str:
    '''
    minor func
    '''
    return 'Supplementary_'+mmcif_file+'.html'


def get_output_file_temp_html(mmcif_file: str) -> str:
    '''
    minor func
    '''
    return 'temp.html'


def get_output_file_pdf(mmcif_file: str) -> str:
    '''
    minor func
    '''
    return mmcif_file+'.pdf'


def get_output_file_json(mmcif_file: str) -> str:
    '''
    minor func
    '''
    return 'ValidationReport_'+mmcif_file+'.json'


def get_supp_file_pdf(mmcif_file: str) -> str:
    '''
    minor func
    '''
    return 'Supplementary_'+mmcif_file+'.pdf'


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
    # print (data_dict)
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


def exv_readable_format(exv: dict) -> list:
    '''
    format exv for supplementary/summary table
    '''

    fin_string = []
    print(exv)
    for ind, el in enumerate(exv['Models']):
        fin_string.append('Model-'+str(el)+': '+'Number of violations = ' +
                          str(exv['Number of violations'][ind]) + ' ')
    return fin_string


def mp_readable_format(mp: dict) -> list:
    '''
    format molprobity resukts for supplementary/summary table
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
    format cx-ms data for supplementary/summary table
    '''

    fin_cx = []
    count = 0
    for key, val in cx_dict.items():
        count += 1
        fin_cx.append('CX-MS Fit of medioid: model # ' + str(count) +
                      ', percentage satisfied ' + str(round(val, 2))+'%')
    return fin_cx


def clean_all():
    '''
    delete all generated files
    '''

    # dirname_ed = os.getcwd()
    os.listdir('.')
    for item in os.listdir('.'):
        if item.endswith('.txt'):
            os.remove(item)
        if item.endswith('.csv'):
            os.remove(item)
        if item.endswith('.json'):
            os.remove(item)
        if item.endswith('.sascif'):
            os.remove(item)


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
        raise (f'Wrong value: {value}. '
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

    if lower > 0:
        oom = order_of_magnitude(lower)
        # Do not allow the range to go below zero
        lower = max(0, lower * (1 - 10 ** (-oom)))

    assert lower >= 0 and upper > 0

    return (lower, upper)
