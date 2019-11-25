###################################
# Script to run validation and write
# to HTML and PDF 
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
from get_input_information import *
from get_molprobity import *
from get_excluded_volume import *
from get_plots import *
import pdfkit
import datetime,time
from flask import Flask, render_template
from multiprocessing import Process, Queue, Pool, Manager
#############################
# Features to add:
# 1) Flask and server support
#############################

config = pdfkit.configuration(wkhtmltopdf='/home/ganesans/PDB-dev/master/pyext/wkhtmltox/bin/wkhtmltopdf')
options = {
    'page-size': 'A4',
    'margin-top': '0.75in',
    'margin-right': '0.75in',
    'margin-bottom': '0.75in',
    'margin-left': '0.75in',
    'enable-javascript': None,
    'javascript-delay':'500'
}
sys.path.append('/home/ganesans/PDB-dev/master/pyext/src/table')
sys.path.append('/home/ganesans/PDB-dev/master/pyext/src/table/images')
d=datetime.datetime.now();t=pytz.timezone("America/Los_Angeles");d1=t.localize(d)
timestamp=d1.strftime("%B %d, %Y --  %I:%M %p")

# Create directory
dirName = 'Output'
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
template_file = "template.html"
#template = templateEnv.get_template(TEMPLATE_FILE)
template_file2 = "template_pdf.html"
#template_pdf = templateEnv.get_template(TEMPLATE_FILE2)
Template_Dict={}
#################################################
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

def get_output_file_temp_html(mmcif_file):
    return 'temp.html'

def get_output_file_pdf(mmcif_file):
    return 'ValidationReport_'+mmcif_file+'.pdf'


all_files=get_all_files(path_dir)

#################################################
# Run one file
################################################

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
    Template_Dict['ID']=I.get_id()
    Template_Dict['Molecule']=I.get_struc_title()
    Template_Dict['Title']=I.get_title()
    Template_Dict['Authors']=I.get_authors()
    Template_Dict['Entry_list']=dict_to_JSlist(I.get_composition())
    Template_Dict['number_of_molecules']=I.get_number_of_models()
    Template_Dict['model_names']=I.get_model_names()
    Template_Dict['number_of_software']=I.get_software_length()
    Template_Dict['soft_list']=dict_to_JSlist(I.get_software_comp())
    Template_Dict['number_of_datasets']=I.get_dataset_length()
    Template_Dict['Datasets_list']=dict_to_JSlist(I.get_dataset_comp())
    Template_Dict['Protocols_number']=I.get_protocol_number()
    Template_Dict['Sampling_list']=dict_to_JSlist(I.get_sampling())
    #print (Template_Dict)
    #d['entry']=Template_Dict
    #return Template_Dict
    
def run_model_quality(mmcif_file):
    I = get_input_information(mmcif_file)
    if I.check_sphere()<1:
        I_mp=get_molprobity_information(mmcif_file)
        if I_mp.check_for_molprobity():
            d_mp=manager.dict()
            runInParallel(I_mp.run_clashscore(d_mp),I_mp.run_ramalyze(d_mp),I_mp.run_rotalyze(d_mp),I_mp.run_molprobity(d_mp))
            a,b=I_mp.process_molprobity(d_mp['molprobity'])
            Template_Dict['bond']=len(a); Template_Dict['angle']=len(b)
            d,e,f=I_mp.get_data_for_quality_at_glance(d_mp['molprobity'])
            Template_Dict['molp_b']=dict_to_JSlist(I_mp.molprobity_detailed_table_bonds(a))
            Template_Dict['molp_a']=dict_to_JSlist(I_mp.molprobity_detailed_table_angles(b))
            Template_Dict['rotascore']=dict_to_JSlist(I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
            Template_Dict['rotalist']=dict_to_JSlist(I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota'])))
            Template_Dict['ramascore']=dict_to_JSlist(I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
            Template_Dict['ramalist']=dict_to_JSlist(I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama'])))
            clashscores,Template_Dict['tot']=I_mp.clash_summary_table(d_mp['clash'])
            Template_Dict['clashscore_list']=dict_to_JSlist(clashscores)
            Template_Dict['clashlist']=I_mp.clash_detailed_table(d_mp['clash'])
            I_mp_plt=plots_mp(mmcif_file)
            I_mp_plt.plot_quality_at_glance_mp(d,e,f)
    else:
        I_ev=get_excluded_volume(mmcif_file)
        model_dict=I_ev.get_all_spheres()
        I_ev.get_exc_vol_for_models(model_dict)
        exv_data=I_ev.get_exc_vol_for_models(I_ev.get_all_spheres())
        Template_Dict['excluded_volume']=dict_to_JSlist(exv_data)
        print (exv_data)
        I_exv_plt=plots_exv(mmcif_file)
        I_exv_plt.plot_quality_at_glance_ev(exv_data)

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

def write_pdf_alone(Template_Dict, template_file):
   pass 

def write_pdf_from_html(html_file):
    pass

def clean_all():
    dirname_ed='/home/ganesans/PDB-dev/master/pyext/src/validation/'    
    os.listdir(dirname_ed)
    for item in os.listdir(dirname_ed):
        if item.endswith('.txt'):
            os.remove(item)


if __name__ == "__main__":
    clean_all()
    manager = Manager() # create only 1 mgr
    d = manager.dict() # create only 1 dict
    start=time.process_time()
    args=[mmcif_file for mmcif_file in all_files]
    run_entry_composition(sys.argv[1])
    run_model_quality(sys.argv[1])
    write_html(sys.argv[1],Template_Dict,template_file)
    #write_pdf(sys.argv[1],Template_Dict,template_file)
    #for filename in all_files:
    #    mmcif_file=filename
    #    mmcif_output = os.path.basename(filename)
    #    run_entry_composition(filename)
    #    run_model_quality(filename)
    #    write_html(mmcif_output,Template_Dict,template_file)
    #    write_pdf(mmcif_output,Template_Dict,template_file)
       # clean_all()
