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
import pdfkit
import datetime,time
import pickle
from flask import Flask, render_template
from multiprocessing import Process, Queue, Pool, Manager
from collections import Counter
import argparse
import json
#############################
# Features to add:
# 1) Flask and server support
#############################

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
config = pdfkit.configuration(wkhtmltopdf='/home/ganesans/PDB-dev/master/pyext/wkhtmltox/bin/wkhtmltopdf')
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
print (dirName)
try:
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")

dirName_supp = 'Supplementary'
try:
    os.mkdir(dirName_supp)
    print("Directory " , dirName_supp ,  " Created ")
except FileExistsError:
    print("Directory " , dirName_supp ,  " already exists")

path_dir='/home/ganesans/PDB-dev/master/pyext/src/data/*'
templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
template_html_main = "main.html"
template_html_data = "data_quality.html"
template_html_model = "model_quality.html"
template_html_comp = "model_composition.html"
template_html_for = "formodeling.html"
template_html_not = "notformodeling.html"
template_html_un = "uncertainty.html"
template_html = "template.html"
template_pdf = "template_pdf.html"
template_file_supp= "supplementary_template.html"
Template_Dict={}

#########################################################################################################################################
# basic functions
#######################################################################################
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

all_files=get_all_files(path_dir)

##########################################################################################################################
# Run one file
################################################

def run_entry_composition(mmcif_file):
    start=time.process_time()
    name=mmcif_file.split('.')[0].split('_')[0]
    I = validation.get_input_information(mmcif_file)
    if I.get_ensembles():
        ensemble_info=validation.dict_to_JSlist(I.get_ensembles())
    else:
        ensemble_info=None
    Template_Dict['ensemble_info']=ensemble_info
    Template_Dict['sphere']=I.check_sphere()
    Template_Dict['num_ensembles']=I.check_ensembles()
    RB,flex,RB_nos,all_nos=I.get_RB_flex_dict()
    Template_Dict['Rigid_Body']=RB_nos
    Template_Dict['Flexible_Unit']=all_nos-RB_nos
    Template_Dict['RB_list']=validation.dict_to_JSlist_rows(RB,flex)
    Template_Dict['RB']=get_RB(validation.dict_to_JSlist_rows(RB,flex))
    Template_Dict['flex']=get_flex(validation.dict_to_JSlist_rows(RB,flex))
    Template_Dict['date']=timestamp
    Template_Dict['ID']=I.get_id()
    Template_Dict['ID_w']=I.get_id().split()
    Template_Dict['ID_T']=I.get_id()[0:6]+'_'+I.get_id()[6:]
    Template_Dict['ID_R']=(I.get_id()[0:6]+'_'+I.get_id()[6:]).split()
    Template_Dict['Molecule']=I.get_struc_title()
    Template_Dict['Title']=I.get_title()
    Template_Dict['Authors']=I.get_authors()
    Template_Dict['Entry_list']=validation.dict_to_JSlist(I.get_composition())
    Template_Dict['number_of_molecules']=I.get_number_of_models()
    Template_Dict['model_names']=I.get_model_names()
    Template_Dict['number_of_software']=I.get_software_length()
    Template_Dict['soft_list']=validation.dict_to_JSlist(I.get_software_comp())
    Template_Dict['number_of_datasets']=I.get_dataset_length()
    Template_Dict['Datasets_list']=validation.dict_to_JSlist(I.get_dataset_comp())
    Template_Dict['Protocols_number']=I.get_protocol_number()
    Template_Dict['Sampling_list']=validation.dict_to_JSlist(I.get_sampling())
    Template_Dict['num_chains']=int(len(I.get_composition()['Chain ID']))/int(len(list(Counter(I.get_composition()['Model ID']).keys())))
    return Template_Dict


def run_model_quality(mmcif_file):
    I = validation.get_input_information(mmcif_file)
    print ("exo",I.check_sphere())
    if I.check_sphere()<1:
        I_mp=get_molprobity_information.get_molprobity_information(mmcif_file)
        if I_mp.check_for_molprobity():
            dirname=os.path.dirname(os.path.abspath(__file__))
            filename = os.path.abspath(os.path.join(dirname, 'static/results/',str(Template_Dict['ID'])+'_temp_mp.txt'))
            print (filename)
            #file_mp=os.getcwd()+'/static/results/'+str(Template_Dict['ID'])+'_temp_mp.txt'
            if os.path.exists(filename):
                d_mp={}
                print ("Molprobity analysis file already exists...\n...assuming clashscores, Ramachandran and rotamer outliers have already been calculated")
                with open(filename,'rb') as fp:
                        d_mp['molprobity']=pickle.load(fp)
                f_rota=os.path.abspath(os.path.join(dirname, 'static/results/',str(Template_Dict['ID'])+'_temp_rota.txt'))
                with open(f_rota,'rb') as fp:
                        d_mp['rota']=pickle.load(fp)
                f_rama=os.path.abspath(os.path.join(dirname, 'static/results/',str(Template_Dict['ID'])+'_temp_rama.txt'))
                with open(f_rama,'rb') as fp:
                        d_mp['rama']=pickle.load(fp)
                f_clash=os.path.abspath(os.path.join(dirname, 'static/results/',str(Template_Dict['ID'])+'_temp_clash.txt'))
                with open(f_clash,'rb') as fp:
                        d_mp['clash']=pickle.load(fp)
            else:
                print ("Molprobity analysis is being calculated...")
                manager = Manager()
                d_mp=manager.dict()
                runInParallel(I_mp.run_clashscore(d_mp),I_mp.run_ramalyze(d_mp),I_mp.run_rotalyze(d_mp),I_mp.run_molprobity(d_mp))
            a,b=I_mp.process_molprobity(d_mp['molprobity'])
            Template_Dict['bond']=len(a); Template_Dict['angle']=len(b)
            d,e,f=I_mp.get_data_for_quality_at_glance(d_mp['molprobity'])
            Template_Dict['molp_b']=validation.dict_to_JSlist(I_mp.molprobity_detailed_table_bonds(a))
            Template_Dict['molp_a']=validation.dict_to_JSlist(I_mp.molprobity_detailed_table_angles(b))
            Template_Dict['rotascore']=validation.dict_to_JSlist(I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
            Template_Dict['rotalist']=validation.dict_to_JSlist(I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota'])))
            Template_Dict['ramascore']=validation.dict_to_JSlist(I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
            Template_Dict['ramalist']=validation.dict_to_JSlist(I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama'])))
            clashscores,Template_Dict['tot']=I_mp.clash_summary_table(d_mp['clash'])
            Template_Dict['clashscore_list']=validation.dict_to_JSlist(clashscores)
            Template_Dict['clashlist']=I_mp.clash_detailed_table(d_mp['clash'])
            I_mp_plt=validation.get_plots.plots_mp(mmcif_file)
            I_mp_plt.plot_quality_at_glance_mp(d,e,f)
            Template_Dict['assess_atomic_segments']='Clashscore: '+ str(d) + ', Ramachandran outliers: '+ str(e)+ ', Sidechain outliers: '+str(f)
            Template_Dict['assess_excluded_volume']=['Not applicable']

    else:
        Template_Dict['assess_atomic_segments']='Not applicable'
        file=os.getcwd()+'/static/results/'+str(Template_Dict['ID'])+'exv.txt'
        if os.path.exists(file):
            print ("Excluded volume file already exists...")
            with open(file, 'r+') as inf:
                line=[ln.replace('[','').replace(']','').replace(',','').split() for ln in inf.readlines()]
            exv_data={'Models':line[0],'Excluded Volume Satisfaction':line[1], 'Number of violations':line[2]}
        else:    
            print ("Excluded volume is being calculated...")
            I_ev=get_excluded_volume.get_excluded_volume(mmcif_file)
            model_dict=I_ev.get_all_spheres()
            exv_data=I_ev.run_exc_vol_parallel(model_dict)
        
        Template_Dict['excluded_volume']=validation.dict_to_JSlist(exv_data)
        I_exv_plt=validation.get_plots.plots_exv(mmcif_file)
        I_exv_plt.plot_quality_at_glance_ev(exv_data)
        Template_Dict['assess_excluded_volume']=exv_readable_format(exv_data)
        print (Template_Dict['assess_excluded_volume'])

def run_sas_validation(mmcif_file):
    I = validation.get_input_information(mmcif_file)
    if I.check_for_sas():
        Template_Dict['sas']=["True"]
        I_sas=sas_validation.sas_validation(mmcif_file)
        Template_Dict['sasdb_code']=I_sas.get_SASBDB_code()
        Template_Dict['parameters_volume']=validation.dict_to_JSlist(I_sas.get_parameters_vol())
        Template_Dict['parameters_mw']=validation.dict_to_JSlist(I_sas.get_parameters_mw())
        Template_Dict['pddf_software']=str(I_sas.get_pddf_software())
        Template_Dict['pddf_dmax']=I_sas.get_pddf_dmax()
        Template_Dict['pddf_rg']=I_sas.get_pddf_rg()
        Template_Dict['number_of_fits']=I_sas.get_number_of_fits()
        Template_Dict['chi_table']=validation.dict_to_JSlist(I_sas.get_chi_table())
        Template_Dict['rg_table']=validation.dict_to_JSlist(I_sas.get_rg_table())
        Template_Dict['r2']=I_sas.get_fit_r2()
        I_sas_plt=validation.sas_validation_plots.sas_validation_plots(mmcif_file)
        I_sas_plt.plot_intensities()
        I_sas_plt.plot_intensities_log()
        I_sas_plt.plot_kratky()
        I_sas_plt.plot_porod_debye()
        I_sas_plt.plot_pddf()
        I_sas_plt.Guinier_plot_fit()
        I_sas_plt.Guinier_plot_residuals()
        if Template_Dict['number_of_fits']>0:
            I_sas_plt.plot_fit()
            I_sas_plt.plot_residuals()
            I_sas.get_fit_image()


def run_supplementary_table(args):
    I = validation.get_input_information(args.f)
    if (I.get_ensembles() is not None) and  (all_same(I.get_ensembles()['Clustering method'])):
        Template_Dict['clustering']=I.get_ensembles()['Clustering method'][0]
    elif I.get_ensembles() is not None:
        Template_Dict['clustering']=', '.join(I.get_ensembles()['Clustering method'])
    else:
        Template_Dict['clustering']='Not applicable'
    Template_Dict['location']=args.ls
    Template_Dict['complex_name']=I.get_struc_title()
    Template_Dict['PDB_ID']=I.get_id()
    Template_Dict['Subunits']=get_subunits(I.get_composition())
    Template_Dict['datasets']=get_datasets(I.get_dataset_details()) if I.get_dataset_details() is not None else 'Not provided or used'
    Template_Dict['physics']=physics
    Template_Dict['software']=get_software(I.get_software_comp())+args.ls
    Template_Dict['struc']=I.get_atomic_coverage()
    Template_Dict['method']=get_method_name(I.get_sampling())
    Template_Dict['method_type']=get_method_type(I.get_sampling())
    Template_Dict['method_details']=args.m
    Template_Dict['models']=', '.join(I.get_ensembles()['Number of models']) if I.get_ensembles() is not None else 'Not applicable' 
    Template_Dict['sampling_validation']=args.sv
    Template_Dict['feature']=I.get_ensembles()['Clustering feature'][0] if I.get_ensembles() is not None else 'Not applicable'
    Template_Dict['validation_input']=args.v1
    Template_Dict['cross_validation']=args.v2
    Template_Dict['model_precision']=', '.join([i+'&#8491' for i in I.get_ensembles()['Cluster precision']]) if I.get_ensembles() is not None else 'Model precision can not be calculated with one structure'
    Template_Dict['restraint_info']=get_restraints_info(I.get_restraints()) if I.get_restraints() is not None else 'Not provided or used'
    Template_Dict['Data_quality']=args.dv
    Template_Dict['clustering']=args.c
    Template_Dict['sampling_validation']=args.sv
    Template_Dict['resolution']=args.res

def write_html(mmcif_file,Template_Dict, template_file):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName,template_file)),"w") as fh:
        fh.write(outputText)

def write_pdf(mmcif_file,Template_Dict, template_file):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName,get_output_file_temp_html(mmcif_file))),"w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName,get_output_file_temp_html(mmcif_file))), os.path.join(os.path.join(dirName,get_output_file_pdf(mmcif_file))) ,configuration=config, options=options)
    os.remove(os.path.join(os.path.join(dirName,get_output_file_temp_html(mmcif_file))))

def write_supplementary_table(mmcif_file,Template_Dict,template_file):
    template = templateEnv.get_template(template_file)
    outputText=template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName_supp,get_supp_file_html(mmcif_file))),"w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName_supp,get_supp_file_html(mmcif_file))), os.path.join(os.path.join(dirName_supp,get_supp_file_pdf(mmcif_file))) ,configuration=config, options=options_supp)

def write_json(mmcif_file,Template_Dict):
    j=json.dumps([{'Category': k, 'Itemized_List': v} for k,v in Template_Dict.items()], indent=4)
    with open(os.path.join(os.path.join(dirName,get_output_file_json(mmcif_file))),"w") as fh:
        fh.write(j)
    fh.close()

def write_pdf_alone(Template_Dict, template_file):
   pass 

def write_pdf_from_html(html_file):
    pass

def clean_all():
    dirname_ed='/home/ganesans/PDB-dev/master/pyext/src/write_report/'    
    os.listdir(dirname_ed)
    for item in os.listdir(dirname_ed):
        if item.endswith('.txt'):
            os.remove(item)

############################################################################################################################
# Run script
#################################################

if __name__ == "__main__":
    #clean_all()
    manager = Manager() # create only 1 mgr
    d = manager.dict() # create only 1 dict
    #start=time.process_time()
    #args=[mmcif_file for mmcif_file in all_files]
    #run_supplementary_table(args)
    run_entry_composition(args.f)
    run_model_quality(args.f)
    #run_supplementary_table(args)
    run_sas_validation(args.f)
    write_html(args.f,Template_Dict,template_html_main)
    write_html(args.f,Template_Dict,template_html_comp)
    write_html(args.f,Template_Dict,template_html_data)
    write_html(args.f,Template_Dict,template_html_model)

 #   write_pdf(args.f,Template_Dict,template_pdf)
    #write_json(args.f,Template_Dict)
    #write_supplementary_table(args.f,Template_Dict,template_file_supp)
        
'''
    for filename in all_files:
        mmcif_file=filename
        print ('calculating validation report for %s'%(os.path.basename(filename)))
        mmcif_output = os.path.basename(filename)
        run_entry_composition(filename)
        run_model_quality(filename)
        write_html(mmcif_output,Template_Dict,template_html)
        write_pdf(mmcif_output,Template_Dict,template_pdf)
        clean_all()
'''    
