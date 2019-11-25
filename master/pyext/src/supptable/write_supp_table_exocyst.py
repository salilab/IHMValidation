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

complex_name='Yeast Exocyst Complex'
PDB_ID='PDBDEV_000000XX'
physics='Excluded volume and Sequence connectivity.'
resolution=['Rigid bodies: 1 residue per bead.','Flexible regions: 50 residues per bead.']
struc='51 %'
method_details=['Replica exchange temperature range used: 1.0-5.0',
                'Number of replicas: 8',
                'Number of runs: 100',
                'Number of structures generated: 2,000,000',
                'Movers for rigid bodies: Random translation up to 4 &#8491, rotation up to 1.0 &#8491',
                'Movers for flexible units: Random translation up to 4 &#8491']
models='9668'
model_precision='38 &#8491 (average RMSF of the solution ensemble with respect to the centroid structure)'
sampling_validation=['Sampling precision : 51 &#8491',
                    'p-value and D-statistic of non-parametric Kolmogorov-Smirnov two-sample test : p-value=0.1 (threshold p-value > 0.05), D-statistic=0.1 (threshold D-statistic < 0.3)',
                    'Homogeneity of proportions &#x3A7<sup>2</sup> test : 0.0/0.098 (threshold p-value > 0.05 or Cramer\'s V <0.1) ']
Number_of_clusters= '1'
Population_of_clusters= '99 %'

validation_input=['Cross-links satisfaction: 98 %', '3D-EM satisfaction/average cross-correlation between model and map densities: 0.82',
                  'Satisfaction of physical principles: 99 % sequence connectivity and 98 % excluded volume' ]

cross_validation= ['3D-EM satisfaction with previously published cryo-EM map (EMDB:6827)/average cross-correlation between model and map densities: 0.78',
                    'Cross-links satisfaction with previously published data (PMID: 29335562): 98 %',
                    '<i>In vivo</i> inter subunit distance satisfaction (PMID: 28129539): 100 %']
location = ['Scripts and data: <a href="https://salilab.org/exocyst">https://salilab.org/exocyst</a>','Output files: <a href="https://zenodo.org/record/XXX">https://zenodo.org/record/XXX</a']             

clustering='distance threshold-based clustering'
ensembles='1'
#################################################
# Run one file
################################################

def run_one(mmcif):
    mmcif_file=mmcif
    name=mmcif_file.split('.')[0].split('_')[0]
    I = get_input_information(mmcif_file)
    RB,flex,RB_nos,all_nos=I.get_RB_flex_dict()
    outputText = template.render(complex_name=complex_name,
                                PDB_ID=PDB_ID,
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
                                struc=struc,
                                method=get_method_name(I.get_sampling()),
                                method_type=get_method_type(I.get_sampling()),
                                method_details=method_details,
                                models=models,
                                ensembles=ensembles,
                                sampling_validation=sampling_validation,
                                Number_of_clusters= Number_of_clusters,
                                clustering=clustering,
                                Population_of_clusters=Population_of_clusters,
                                validation_input=validation_input,
                                cross_validation=cross_validation,
                                model_precision=model_precision,
                                restraint_info=get_restraints_info(I.get_restraints())
                                )
                                            

    with open('Supp_table/Table_'+mmcif_file+'.html',"w") as fh:
        fh.write(outputText)
        
    pdfkit.from_file('Supp_table/Table_'+mmcif_file+'.html', 'Supp_table/Table_'+mmcif_file+'.pdf' ,configuration=config, options=options)



if __name__ == "__main__":
    run_one(sys.argv[1])
