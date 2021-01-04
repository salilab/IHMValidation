###################################
# Script : 
# 1) Contains trivial and useful
# functions
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

import pandas as pd
import sys,os,glob
import numpy as np
import pickle
from collections import Counter
import argparse


def dict_to_JSlist(d):
    '''
    convert dictionary to list of lists
    '''
    L = []
    if bool(d):
        if len(list(d.keys()))>0:
            L.append(list(d.keys()))
            target=list(d.values())
            for i in range(len(target[0])):
                ltt=[]
                for j in target:
                    ltt.append(str(j[i]))
                L.append(ltt)
    return L

def format_RB_text(tex):
    '''
    convert RB information to text for supp table
    '''
    val=''
    for a in tex:
        for b in a:
            if b==a[-1] and a==tex[-1]:
                val+=str(b)+'. '
            elif b==a[-1] and a!= tex[-1]:
                val+=str(b)+', '
            else:
                val+=str(b)+':'
    if val=='':
        val='-'
    return val

def format_flex_text(tex):
    '''
    convert flex information to text for supp table
    '''
    val=''
    for a in tex:
        for b in a:
            if b==a[-1] and a==tex[-1]:
                val+=str(b)+'. '
            else:
                val+=str(b)+', '

    if val=='':
        val='-'

    return val


def format_tupple(tex):
    return str(tex[0])+'-'+str(tex[1])


def dict_to_JSlist_rows(d1,d2):
    '''
    '''
    L=[]
    L.append(['Chain ID','Rigid bodies','Non-rigid segments'])
    for i,j in d1.items():
        L.append([i,format_RB_text(j),format_flex_text(d2[i])])
    return L

def islistempty(inlist):
    if isinstance (inlist,list):
        return all(map(islistempty,inlist))
    return False

def cat_list_string(listn):
    result=' '
    for i in range(len(listn)):
        if i==0:
            result += str(listn[i])
        else:
            result += ','
            result += str(listn[i])
    return result

def get_key_from_val(dict1,val1):
    return dict1.keys()[dict1.values().index(val1)]

def get_val_from_key(dict1,key1):
    return dict1[key1]

def get_name(name):
    #if str(name) in ['?','',1,'.']:
    #    return 'None Listed'
    #else:
    return str(name)

def get_copy(name):
    if str(name)=='?':
        copy='None listed'
    elif '.' in name:
        copy=(name.split('.')[1]).split('@')[0]
    else:
        copy=1
    return copy

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
    return mmcif_file+'.pdf'

def get_output_file_json(mmcif_file):
    return 'ValidationReport_'+mmcif_file+'.json'

def get_supp_file_pdf(mmcif_file):
    return 'Supplementary_'+mmcif_file+'.pdf'

def get_all_files(path_dir):
    return glob.glob(path_dir)

def get_subunits(sub_dict):
    n=len(sub_dict['Model ID'])
    sublist=['%s: Chain %s (%d residues)' % (sub_dict['Subunit name'][i],sub_dict['Chain ID'][i],sub_dict['Total residues'][i]) for i in range(0,n)]
    return list(set(sublist))

def get_datasets(data_dict):
    n=len(data_dict['ID'])
    #print (data_dict)
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
    try:
        dataset=[(restraints['Restraint info'][i],restraints['Restraint type'][i]) for i in range(0,n)]
    except:
        new_restraints=dict()
        for key,val in restraints.items():
            new_restraints[key]=list(set(val))
        n=min(len(new_restraints['Restraint info']),len(new_restraints['Restraint type']))
        dataset=[(new_restraints['Restraint info'][i],new_restraints['Restraint type'][i]) for i in range(0,n)]

    for i,j in Counter(dataset).items():
        datalist.append(['%s unique %s: %s' % (j,i[1],i[0])])
    print (datalist)
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
    #print (exv)
    for i,j in enumerate(exv['Models']):
        fin_string.append('Model-'+str(j)+': '+'Number of violations-' + str(exv['Number of violations'][i]) + ' ')
    return fin_string

def get_rg_data(rg_dict):
    fin_rg=[]
    for key,val in rg_dict.items():
        fin_rg.append(key+': Rg from Gunier is '+ str(val[0])+'nm and Rg from p(r) is '+ str(val[1])+'nm')
    return fin_rg

def get_rg_data_fits(rg_dict):
    fin_rg=[]
    for key,val in rg_dict.items():
        for i,j in enumerate(val):
            count=i+1
            fin_rg.append(key+': Fit '+ str(count) +' with &#x3A7;&#xb2; value '+ str(j))
    return fin_rg

def get_cx_data_fits(cx_dict):
    fin_cx=[];count=0
    for key,val in cx_dict.items():
        count+=1
        fin_cx.append('CX-MS Fit of medioid: model # '+ str(count)+', percentage satisfied ' +str(round(val,2))+'%')
    return fin_cx


def clean_all():
    #dirname_ed='/Users/saijananiganesan/Desktop/PDB-dev/working'
    #dirname_ed=os.path.normpath(os.getcwd() + os.sep + os.pardir)
    dirname_ed=os.getcwd()
    #print ("dirname",dirname_ed)
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


