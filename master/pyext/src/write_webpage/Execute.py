###################################
# Script to run validation and write
# to HTML and PDF 
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import pytz
import jinja2
import pandas as pd
import sys,os,glob
import numpy as np
import validation
from validation import get_excluded_volume
from validation import get_molprobity_information
from validation import get_plots,sas_validation,sas_validation_plots
#import pdfkit
import datetime,time
import pickle
from multiprocessing import Process, Queue, Pool, Manager
from collections import Counter
import argparse
import json
#import validation.utility
from validation.WriteReport import WriteReport

####################################################################################################################
# Add input arguments for supp table
#####################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=str, default='Yes', help ="Physical principles used in modeling yes/no?")
parser.add_argument('-f',default='PDBDEV_00000001.cif',help ="Input mmcif file")
parser.add_argument('-ls', type=list, default=['No location specified'], help ="add location of your scripts")
parser.add_argument('-ld', type=list, default=['No location specified'], help ="add location of your analysis files")
parser.add_argument('-m', type=list, default=['Method details unspecified'], help ="add information on your method")
parser.add_argument('-models', type=str, default='1', help ="number of models in an ensemble, if you have multiple ensembles, add comma-separated string")
parser.add_argument('-c', type=str, default='Distance threshold-based clustering', help ="The type of clustering algorithm used to analyze the ensemble")
parser.add_argument('-mp', type=str, default='10 &#8491 (average RMSF of the solution ensemble with respect to the centroid structure)', help ="add model precision. Model precision is defined as average RMSF of the solution ensemble with respect to the centroid structure")
parser.add_argument('-sv', type=list, default=['Information related to sampling validation has not been provided' ], help ="add model precision. Model precision is defined as average RMSF of the solution ensemble with respect to the centroid structure")
parser.add_argument('-v1', type=list, default=['Fit of model to information used to compute it has not been determined' ], help ="Add information on satisfaction of input data/restraints")
parser.add_argument('-v2', type=list, default=['Fit of model to information not used to compute it has not been determined' ], help ="Add information on satisfaction of data not used for modeling")
parser.add_argument('-dv', type=list, default=['Quality of input data has not be assessed' ], help ="Add information on quality of input data")
parser.add_argument('-res', type=list, default=['Rigid bodies: 1 residue per bead.','Flexible regions: 50 residues per bead.'], help ="Add information on model quality (molprobity or excluded volume)")

args = parser.parse_args()
if args.p.upper() == 'YES':
    physics='Excluded volume and Sequence connectivity.'
else:
    physics='Physical principles were not used for modeling'
#############################################################################################################################
# Input for Jinja
####################################################################################
#config = pdfkit.configuration(wkhtmltopdf='/home/ganesans/PDB-dev/master/pyext/wkhtmltox/bin/wkhtmltopdf')
options = {
    'page-size': 'A4',
    'margin-top': '0.75in',
    'margin-right': '0.75in',
    'margin-bottom': '0.75in',
    'margin-left': '0.75in',
    'enable-javascript': None,
    'javascript-delay':'500',
    'header-left':'[page] of [topage]',
    'footer-center':'Full wwPDB IM Structure Validation Report',
    'footer-line':'',
    'header-line':'',
    'footer-spacing':'5',
    'header-spacing':'5'
}

options_supp = {
    'page-size': 'A4',
    'margin-top': '0.75in',
    'margin-right': '0.75in',
    'margin-bottom': '0.75in',
    'margin-left': '0.75in',
    'enable-javascript': None,
    'javascript-delay':'500',
    'header-left':'[page] of [topage]',
    'footer-center':'wwPDB IM Methods Table',
    'footer-line':'',
    'header-line':'',
    'footer-spacing':'5',
    'header-spacing':'5'
}

sys.path.append('/home/ganesans/PDB-dev/master/pyext/src/table')
sys.path.append('/home/ganesans/PDB-dev/master/pyext/src/table/images')
d=datetime.datetime.now();t=pytz.timezone("America/Los_Angeles");d1=t.localize(d)
timestamp=d1.strftime("%B %d, %Y --  %I:%M %p")

# Create directory
dirName = 'Output/'+str(args.f).split('.')[0]
dirName_flask = 'FlaskApp/templates/'+str(args.f).split('.')[0]

print (dirName)

try:
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")

try:
    os.mkdir(dirName_flask)
    print("Directory " , dirName_flask ,  " Created ") 
except FileExistsError:
    print("Directory " , dirName_flask,  " already exists")

dirName_supp = 'Supplementary'

try:
    os.mkdir(dirName_supp)
    print("Directory " , dirName_supp ,  " Created ")
except FileExistsError:
    print("Directory " , dirName_supp ,  " already exists")

path_dir='/home/ganesans/PDB-dev/master/pyext/src/data/*'
templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
template_html_main = "templates/main.html"
template_html_data = "templates/data_quality.html"
template_html_model = "templates/model_quality.html"
template_html_comp = "templates/model_composition.html"
template_html_for = "templates/formodeling.html"
template_html_not = "templates/notformodeling.html"
template_html_un = "templates/uncertainty.html"

template_flask_main = "templates_flask/main.html"
template_flask_data = "templates_flask/data_quality.html"
template_flask_model = "templates_flask/model_quality.html"
template_flask_comp = "templates_flask/model_composition.html"
template_flask_for = "templates_flask/formodeling.html"
template_flask_not = "templates_flask/notformodeling.html"
template_flask_un = "templates_flask/uncertainty.html"

template_html = "template.html"
template_pdf = "template_pdf.html"
template_file_supp= "supplementary_template.html"
Template_Dict={}
Template_Dict['date']=timestamp

all_files=validation.utility.get_all_files(path_dir)

#############################################################################################################################
# Jinja scripts
#############################################################################################################################

def write_html(Template_Dict, template_file,dirName):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    template_file=template_file.split('/')[1]
    with open(os.path.join(os.path.join(dirName,template_file)),"w") as fh:
        fh.write(outputText)

def write_pdf(mmcif_file,Template_Dict, template_file,dirName):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName,validation.utility.get_output_file_temp_html(mmcif_file))),"w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName,validation.utility.get_output_file_temp_html(mmcif_file))), 
        os.path.join(os.path.join(dirName,validation.utility.get_output_file_pdf(mmcif_file))) ,
        configuration=config, 
        options=options)
    os.remove(os.path.join(os.path.join(dirName,validation.utility.get_output_file_temp_html(mmcif_file))))

def write_supplementary_table(mmcif_file,Template_Dict,template_file,dirName_supp):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName_supp,validation.utility.get_supp_file_html(mmcif_file))),"w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName_supp,validation.utility.get_supp_file_html(mmcif_file))), 
        os.path.join(os.path.join(dirName_supp,validation.utility.get_supp_file_pdf(mmcif_file))) ,
        configuration=config, 
        options=options_supp)

def write_json(mmcif_file,Template_Dict):
    j=json.dumps([{'Category': k, 'Itemized_List': v} for k,v in Template_Dict.items()], indent=4)
    with open(os.path.join(os.path.join(dirName,validation.utility.get_output_file_json(mmcif_file))),"w") as fh:
        fh.write(j)
    fh.close()
############################################################################################################################
# Run script
#################################################

if __name__ == "__main__":
    validation.utility.clean_all()

    manager = Manager() # create only 1 mgr
    d = manager.dict() # create only 1 dict
    report=WriteReport(args.f)
    template_dict=report.run_entry_composition(Template_Dict)
    print ("data",template_dict['Data'])
    template_dict,clashscore,rama,sidechain,exv_data=report.run_model_quality(template_dict)
    template_dict,sas_data,sas_fit=report.run_sas_validation(template_dict)
    print (clashscore,rama,sidechain,exv_data)
    print ("sasfit",sas_fit)
    report.run_quality_glance(clashscore,rama,sidechain,exv_data,sas_data,sas_fit)

    write_html(template_dict,template_html_main,dirName)
    write_html(template_dict,template_html_comp,dirName)
    write_html(template_dict,template_html_data,dirName)
    write_html(template_dict,template_html_model,dirName)
    write_html(template_dict,template_html_for,dirName)
    write_html(template_dict,template_html_not,dirName)
    write_html(template_dict,template_html_un,dirName)

    write_html(template_dict,template_flask_main,dirName_flask)
    write_html(template_dict,template_flask_comp,dirName_flask)
    write_html(template_dict,template_flask_data,dirName_flask)
    write_html(template_dict,template_flask_model,dirName_flask)
    write_html(template_dict,template_flask_for,dirName_flask)
    write_html(template_dict,template_flask_not,dirName_flask)
    write_html(template_dict,template_flask_un,dirName_flask)


    validation.utility.clean_all()




