import pandas as pd
import sys,os,glob
import numpy as np
import pickle
from collections import Counter
import argparse

def get_all_files(path_dir):
    return glob.glob(path_dir)

def runInParallel(*fns):
  proc = []
  for fn in fns:
    p = Process(target=fn,args=d)
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

def get_output_file_html(mmcif_file):
    return 'ValidationReport_'+mmcif_file+'.html'

def get_supp_file_html(mmcif_file):
    return 'Supplementary_'+mmcif_file+'.html'

def get_output_file_temp_html(mmcif_file):
    return 'temp.html'

def get_output_file_pdf(mmcif_file):
    return 'ValidationReport_'+mmcif_file+'.pdf'

def get_output_file_json(mmcif_file):
    return 'ValidationReport_'+mmcif_file+'.json'

def get_supp_file_pdf(mmcif_file):
    return 'Supplementary_'+mmcif_file+'.pdf'

def get_all_files(path_dir):
    return glob.glob(path_dir)

def get_subunits(sub_dict):
    n=len(sub_dict['Model ID'])
    sublist=['%s: Chain %s (%d residues)' % (sub_dict['Subunit name'][i],sub_dict['Chain ID'][i],sub_dict['Total residues'][i]) for i in range(0,n)]
    return sublist

def get_datasets(data_dict):
    n=len(data_dict['ID'])
    print (data_dict)
    datalist=['%s, %s' % (data_dict['Dataset type'][i],data_dict['Details'][i]) for i in range(0,n)]
    return datalist

def get_software(data_dict):
    if len(data_dict)>0:
        n=len(data_dict['ID'])
        datalist=['%s (version %s)' % (data_dict['Software name'][i],data_dict['Software version'][i]) for i in range(0,n)]
        return datalist
    else:
        return ['Software details not provided']

def get_RB(data_list):
    n=len(data_list)
    datalist=['%s: %s ' % (data_list[i][0],data_list[i][1],) for i in range(1,n)]
    return datalist

def get_flex(data_list):
    n=len(data_list)
    datalist=['%s: %s ' % (data_list[i][0],data_list[i][2],) for i in range(1,n)]
    return datalist

def get_method_name(sample_dict):
    datalist='%s ' % (sample_dict['Method name'][0])
    return datalist.replace('monte carlo','Monte Carlo')

def get_method_type(sample_dict):
    datalist='%s ' % (sample_dict['Method type'][0])
    return datalist.replace('monte carlo','Monte Carlo')

def get_restraints_info(restraints):
    n=len(restraints['Restraint type'])
    datalist=[]
    dataset=[(restraints['Restraint info'][i],restraints['Restraint type'][i]) for i in range(0,n)]
    for i,j in Counter(dataset).items():
        datalist.append(['%s unique %s: %s' % (j,i[1],i[0])])
    return datalist

def format_list_text(sublist):
    val=''
    for a in sublist:
        if a==sublist[-1]:
            val+=str(a)+'. '
        else:
            val+=str(a)+', '
    if val =='':
        val='-'
    return val

def all_same(items):
    return all(x==items[0] for x in items)

def exv_readable_format(exv):
    fin_string=[]
    print (exv)
    for i,j in enumerate(exv['Models']):
        fin_string.append('Model-'+str(j)+': '+'Number of violations-' + str(exv['Number of violations'][i]) + ' ')
    return fin_string

def clean_all():
    #dirname_ed='/Users/saijananiganesan/Desktop/PDB-dev/working'
    #dirname_ed=os.path.normpath(os.getcwd() + os.sep + os.pardir)
    dirname_ed=os.getcwd()
    print ("dirname",dirname_ed)
    os.listdir(dirname_ed)
    for item in os.listdir(dirname_ed):
        if item.endswith('.txt'):
            os.remove(item)
        if item.endswith('.csv'):
            os.remove(item)
        if item.endswith('.json'):
            os.remove(item)
        if item.endswith('.sascif'):
            os.remove(item)


