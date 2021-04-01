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


def dict_to_JSlist(d: dict) -> list:
    '''
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
                sublist.append(str(el[ind]))
            output_list.append(sublist)
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

    if val == '':
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

    if val == '':
        val = ['-']

    return ''.join(val)


def format_tuple(tex: list) -> str:
    return str(tex[0])+'-'+str(tex[1])


def dict_to_JSlist_rows(dict1: dict, dict2: dict) -> list:
    output_list = []
    output_list.append(['Chain ID', 'Rigid bodies', 'Non-rigid segments'])
    for ind, el in dict1.items():
        output_list.append(
            [ind, format_RB_text(el), format_flex_text(dict2[ind])])
    return output_list


def islistempty(inlist: list) -> bool:
    if isinstance(inlist, list):
        return all(map(islistempty, inlist))
    return False


def cat_list_string(listn: list) -> str:
    result = ' '
    for ind in range(len(listn)):
        if ind == 0:
            result += str(listn[ind])
        else:
            result += ','
            result += str(listn[ind])
    return result


def get_key_from_val(dict1: dict, val1: str) -> list:
    return dict1.keys()[dict1.values().index(val1)]


def get_val_from_key(dict1: dict, key1: str) -> list:
    return dict1[key1]


def get_name(name):
    return str(name)


def get_copy(name):
    if str(name) == '?':
        copy = 'None listed'
    elif '.' in name:
        copy = (name.split('.')[1]).split('@')[0]
    else:
        copy = 1
    return copy


def get_all_files(path_dir):
    return glob.glob(path_dir)


def runInParallel(*fns, d):
    proc = []
    for fn in fns:
        p = Process(target=fn, args=d)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()


def runInParallel_noargs(*fns):
    proc = []
    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()


def get_output_file_html(mmcif_file: str) -> str:
    return 'ValidationReport_'+mmcif_file+'.html'


def get_supp_file_html(mmcif_file: str) -> str:
    return 'Supplementary_'+mmcif_file+'.html'


def get_output_file_temp_html(mmcif_file: str) -> str:
    return 'temp.html'


def get_output_file_pdf(mmcif_file: str) -> str:
    return mmcif_file+'.pdf'


def get_output_file_json(mmcif_file: str) -> str:
    return 'ValidationReport_'+mmcif_file+'.json'


def get_supp_file_pdf(mmcif_file: str) -> str:
    return 'Supplementary_'+mmcif_file+'.pdf'


def get_subunits(sub_dict: dict) -> list:
    model_number = len(sub_dict['Model ID'])
    sublist = ['%s: Chain %s (%d residues)' % (sub_dict['Subunit name'][i], sub_dict['Chain ID']
                                               [i], sub_dict['Total residues'][i]) for i in range(model_number)]
    return list(set(sublist))


def get_datasets(data_dict: dict) -> list:
    dataset_number = len(data_dict['ID'])
    # print (data_dict)
    datalist = ['%s, %s' % (data_dict['Dataset type'][i], data_dict['Details'][i])
                for i in range(dataset_number)]
    return datalist


def get_software(data_dict: dict) -> list:
    if len(data_dict) > 0:
        dataset_number = len(data_dict['ID'])
        datalist = ['%s (version %s)' % (data_dict['Software name'][i],
                                         data_dict['Software version'][i]) for i in range(dataset_number)]
        return datalist
    else:
        return ['Software details not provided']


def get_RB(data_list: list) -> list:
    data_num = len(data_list)
    #datalist = ['%s: %s ' % (data_list[i][0], data_list[i][1],)
    #           for i in range(1, data_num)]
    datalist=[]
    for i in range(1, data_num):
        if len(data_list[i][1])<1:
            data_list[i][1]='None'
        datalist.append('%s: %s ' % (data_list[i][0], data_list[i][1]))
    return datalist


def get_flex(data_list: list) -> list:
    data_num = len(data_list)
    datalist = ['%s: %s ' % (data_list[i][0], data_list[i][2],)
                for i in range(1, data_num)]
    return datalist


def get_method_name(sample_dict: dict) -> str:
    datastr = '%s ' % (sample_dict['Method name'][0])
    return datastr.replace('monte carlo', 'Monte Carlo')


def get_method_type(sample_dict: dict) -> str:
    datastr = '%s ' % (sample_dict['Method type'][0])
    return datastr.replace('monte carlo', 'Monte Carlo')


def get_restraints_info(restraints: dict) -> list:
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
    return all(x == items[0] for x in items)


def exv_readable_format(exv: dict) -> list:
    fin_string = []
    for ind, el in enumerate(exv['Models']):
        fin_string.append('Model-'+str(el)+': '+'Number of violations-' +
                          str(exv['Number of violations'][ind]) + ' ')
    return fin_string


def get_rg_data(rg_dict: dict) -> list:
    fin_rg = []
    for key, val in rg_dict.items():
        fin_rg.append(key+': Rg from Gunier is ' +
                      str(val[0])+'nm and Rg from p(r) is ' + str(val[1])+'nm')
    return fin_rg


def get_rg_data_fits(rg_dict: dict) -> list:
    fin_rg = []
    for key, val in rg_dict.items():
        for ind, el in enumerate(val):
            count = ind+1
            fin_rg.append(key+': Fit ' + str(count) +
                          ' with &#x3A7;&#xb2; value ' + str(el))
    return fin_rg


def get_cx_data_fits(cx_dict: dict) -> list:
    fin_cx = []
    count = 0
    for key, val in cx_dict.items():
        count += 1
        fin_cx.append('CX-MS Fit of medioid: model # ' + str(count) +
                      ', percentage satisfied ' + str(round(val, 2))+'%')
    return fin_cx


def clean_all():
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
