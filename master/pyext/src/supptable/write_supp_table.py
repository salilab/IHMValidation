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
import validation
import pdfkit
import datetime
import argparse

################################################
parser = argparse.ArgumentParser()
   
parser.add_argument('-p', type=str, default='Yes', help ="Physical principles used in modeling yes/no?")
parser.add_argument('-f', type=str ,default='PDBDEV_00000001.cif',help ="Input mmcif file")
parser.add_argument('-ls', type=list, default=['No location specified'], help ="add location of your scripts")
parser.add_argument('-ld', type=list, default=['No location specified'], help ="add location of your analysis files")
parser.add_argument('-m', type=list, default=['Method details unspecified'], help ="add information on your method")
parser.add_argument('-sp', type=list, default=['Method details unspecified'], help ="add information on your method")
parser.add_argument('-models', type=str, default='1', help ="number of models in an ensemble, if you have multiple ensembles, add comma-separated string")
parser.add_argument('-mp', type=str, default='10 &#8491 (average RMSF of the solution ensemble with respect to the centroid structure)', help ="add model precision. Model precision is defined as average RMSF of the solution ensemble with respect to the centroid structure")
parser.add_argument('-sv', type=list, default=['Information related to sampling validation has not been provided' ], help ="add model precision. Model precision is defined as average RMSF of the solution ensemble with respect to the centroid structure")
parser.add_argument('-v1', type=list, default=['Fit of model to information used to compute it has not been determined' ], help ="Add information on satisfaction of input data/restraints")
parser.add_argument('-v2', type=list, default=['Fit of model to information not used to compute it has not been determined' ], help ="Add information on satisfaction of data not used for modeling")
parser.add_argument('-dv', type=list, default=['Quality of input data has not be assessed' ], help ="Add information on quality of input data")
parser.add_argument('-mv', type=list, default=['Model quality has not been determined' ], help ="Add information on model quality (molprobity or excluded volume)")

args = parser.parse_args()
print (args)

################################################
# Define variables
################################################

path_dir='/home/ganesans/PDB-dev/master/pyext/src/data/*'
templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "supplementary_template.html"
template = templateEnv.get_template(TEMPLATE_FILE)
Temp_Dict={}
################################################
# Variables from user
################################################

if args.p in ['Yes','yes']:
    physics='Excluded volume and Sequence connectivity.'
else:
    physics=''
resolution=['Rigid bodies: 1 residue per bead.','Flexible regions: 10 residues per bead.']
struc=''

method_details=[]
sampling_validation=[]
validation_input=[]
cross_validation= []
location = []

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
