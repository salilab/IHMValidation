###################################
# Script to write supplementaty table
# from mmcif file
# HTML and PDF will be written
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import pytz
import jinja2
import sys
import pandas as pd
import glob
import sys,os
import numpy as np
from create_report import *
import pdfkit
import datetime
import argparse

################################################
parser = argparse.ArgumentParser()
   
parser.add_argument('-p', type=str, default='Yes', help ="Physical principles used in modeling yes/no?")
parser.add_argument('-f', help ="Input mmcif file")

args = parser.parse_args()
print (args)
################################################
# Basic functions
################################################

def get_all_files(path_dir):
    return glob.glob(path_dir)

def get_subunits(sub_dict):
    n=len(sub_dict['Model ID'])
    sublist=['%s: Chain %s (%d residues)' % (sub_dict['Subunit name'][i],sub_dict['Chain ID'][i],sub_dict['Total residues'][i]) for i in range(0,n)]
    return sublist

def get_datasets(data_dict):
    n=len(data_dict['ID'])
    datalist=['%s, %s' % (data_dict['Dataset type'][i],data_dict['Details'][i]) for i in range(0,n)]
    return datalist

def get_software(data_dict):
    n=len(data_dict['ID'])
    datalist=['%s (version %s)' % (data_dict['Software name'][i],data_dict['Software version'][i]) for i in range(0,n)]
    return datalist

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
    if 'DerivedDistanceRestraint' in restraints['Restraint type']:
        dataset=[restraints['Restraint info'][i][1] for i in range(0,n)]
        unique=np.unique(np.array(dataset))
        nos=[];i=0
        for k,l in enumerate(unique):
            nos.append(len([1 for i in range(0,n) if 'DerivedDistanceRestraint' in restraints['Restraint type'] if l in restraints['Restraint info'][i][1]]))
        datalist=['DerivedDistanceRestraint: %s distinct restraints, %s dataset' % (nos[i], unique[i]) for i in range(0, len(unique))]
    else:
        datalist=['%s: %s, %s' % (restraints['Restraint type'][i], restraints['Restraint info'][i][1],restraints['Restraint info'][i][0]) for i in range(0,n)] 
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

################################################
# Template/text file
################################################

config = pdfkit.configuration(wkhtmltopdf='/home/ganesans/PDB-dev/master/pyext/wkhtmltox/bin/wkhtmltopdf')
options = {
    'page-size': 'A4',
    'margin-top': '0.75in',
    'margin-right': '0.75in',
    'margin-bottom': '0.75in',
    'margin-left': '0.75in',
}
sys.path.append('/home/ganesans/PDB-dev/master/pyext/src/table')
sys.path.append('/home/ganesans/PDB-dev/master/pyext/src/table/images')

dirName = 'Supp_table'
try:
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")

Template_Dict={}
################################################
# Define variables
################################################

path_dir='/home/ganesans/PDB-dev/master/pyext/src/data/*'
templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "supplementary_template.html"
template = templateEnv.get_template(TEMPLATE_FILE)

################################################
# Variables from user
################################################

if args.p in ['Yes','yes']:
    Template_Dict['physics']='Excluded volume and Sequence connectivity.'
else:
    Template_Dict['physics']=''
Template_Dict['resolution']=['Rigid bodies: 1 residue per bead.','Flexible regions: 10 residues per bead.']
Template_Dict['struc']=''

Template_Dict['method_details']=[]
Template_Dict['sampling_validation']=[]
Template_Dict['validation_input']=[]
Template_Dict['cross_validation']= []
Template_Dict['location'] = []

#################################################
# Run one file
################################################
###
# edit get_ensembles_when_no_ensembles_exist
def run_one(mmcif):
    mmcif_file=mmcif
    name=mmcif_file.split('.')[0].split('_')[0]
    I = get_input_information(mmcif_file)
    I.dataset_id_type_dic()
    RB,flex,RB_nos,all_nos=I.get_RB_flex_dict()
    print (RB)
    print (RB_nos)
    print (all_nos)
    if (I.get_ensembles() is not None) and  (all_same(I.get_ensembles()['Clustering method'])):
        clustering=I.get_ensembles()['Clustering method'][0]
    elif I.get_ensembles() is not None:
        clustering=', '.join(I.get_ensembles()['Clustering method'])
    else:
        clustering=''

def run_entry_composition(mmcif_file):
    start=time.process_time()
    name=mmcif_file.split('.')[0].split('_')[0]
    I = get_input_information(mmcif_file)
    if I.get_ensembles():
        ensemble_info=dict_to_JSlist(I.get_ensembles())
    else:
        ensemble_info=None
    Template_Dict['ensemble_info']=ensemble_info
    Template_Dict['sphere']=I.check_sphere()
    Template_Dict['num_ensembles']=I.check_ensembles()
    RB,flex,RB_nos,all_nos=I.get_RB_flex_dict()
    Template_Dict['Rigid_Body']=RB_nos
    Template_Dict['Flexible_Unit']=all_nos-RB_nos
    Template_Dict['RB_list']=dict_to_JSlist_rows(RB,flex)
    Template_Dict['date']=timestamp
    Template_Dict['PDB_ID']=I.get_id()
    Template_Dict['Molecule']=I.get_struc_title()
    Template_Dict['Title']=I.get_title()
    Template_Dict['Authors']=I.get_authors()
    Template_Dict['Subunits']=get_subunits(I.get_composition())
    Template_Dict['number_of_molecules']=I.get_number_of_models()
    Template_Dict['model_names']=I.get_model_names()
    Template_Dict['number_of_software']=I.get_software_length()
    Template_Dict['soft_list']=dict_to_JSlist(I.get_software_comp())
    Template_Dict['number_of_datasets']=I.get_dataset_length()
    Template_Dict['datasets']=get_datasets(I.get_dataset_details())
    Template_Dict['Protocols_number']=I.get_protocol_number()
    Template_Dict['Sampling_list']=dict_to_JSlist(I.get_sampling())
    Template_Dict['complex_name']=I.get_struc_title()


    RB,flex,RB_nos,all_nos=I.get_RB_flex_dict()
    outputText = template.render(complex_name=I.get_struc_title(),
                                software=get_software(I.get_software_comp())+location,
                                RB_nos=RB_nos,
                                flex_nos=all_nos-RB_nos,
                                RB=get_RB(dict_to_JSlist_rows(RB,flex)),
                                flex=get_flex(dict_to_JSlist_rows(RB,flex)),
                                resolution=resolution,
                                struc=I.get_atomic_coverage(),
                                method=get_method_name(I.get_sampling()),
                                method_type=get_method_type(I.get_sampling()),
                                method_details=method_details,
                                models=', '.join(I.get_ensembles()['Number of models']),
                                ensembles=len(I.get_ensembles()['Ensemble number']),
                                sampling_validation=sampling_validation,
                                clustering=clustering,
                                feature=I.get_ensembles()['Clustering feature'][0],
                                validation_input=validation_input,
                                cross_validation=cross_validation,
                                model_precision=', '.join([i+'&#8491' for i in I.get_ensembles()['Cluster precision']]),
                                restraint_info=get_restraints_info(I.get_restraints())
                                )

def write_html(mmcif_file,Template_Dict, template_file):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName,get_output_file_html(mmcif_file))),"w") as fh:
        fh.write(outputText)

def write_pdf(mmcif_file,Template_Dict, template_file):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName,get_output_file_temp_html(mmcif_file))),"w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName,get_output_file_temp_html(mmcif_file))), os.path.join(os.path.join(dirName,get_output_file_pdf(mmcif_file))) ,configuration=config, options=options)
    os.remove(os.path.join(os.path.join(dirName,get_output_file_temp_html(mmcif_file))))



    #print (I.get_restraints()['Restraint type'])
    RB,flex,RB_nos,all_nos=I.get_RB_flex_dict()
    outputText = template.render(complex_name=I.get_struc_title(),
                                PDB_ID=I.get_id(),
                                location=location,
                                Subunits=get_subunits(I.get_composition()),
                                datasets=get_datasets(I.get_dataset_details()),
                                physics=physics,
                                software=get_software(I.get_software_comp())+location,
                                RB_nos=RB_nos,
                                flex_nos=all_nos-RB_nos,
                                RB=get_RB(dict_to_JSlist_rows(RB,flex)),
                                flex=get_flex(dict_to_JSlist_rows(RB,flex)),
                                resolution=resolution,
                                struc=I.get_atomic_coverage(),
                                method=get_method_name(I.get_sampling()),
                                method_type=get_method_type(I.get_sampling()),
                                method_details=method_details,
                                models=', '.join(I.get_ensembles()['Number of models']),
                                ensembles=len(I.get_ensembles()['Ensemble number']),
                                sampling_validation=sampling_validation,
                                clustering=clustering,
                                feature=I.get_ensembles()['Clustering feature'][0],
                                validation_input=validation_input,
                                cross_validation=cross_validation,
                                model_precision=', '.join([i+'&#8491' for i in I.get_ensembles()['Cluster precision']]),
                                restraint_info=get_restraints_info(I.get_restraints())
                                )
                                            

    with open('Supp_table/Table_'+mmcif_file+'.html',"w") as fh:
        fh.write(outputText)
        
    pdfkit.from_file('Supp_table/Table_'+mmcif_file+'.html', 'Supp_table/Table_'+mmcif_file+'.pdf' ,configuration=config, options=options)



if __name__ == "__main__":
    run_one(args.f)
