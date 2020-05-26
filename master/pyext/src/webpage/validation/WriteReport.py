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
from validation import utility

class WriteReport(object):
	def __init__(self,mmcif_file):
		self.mmcif_file = mmcif_file
		self.I=validation.get_input_information(self.mmcif_file)

	def run_entry_composition(self,Template_Dict):
		start=time.process_time()
		name=self.mmcif_file.split('.')[0].split('_')[0]
		if self.I.get_ensembles():
			ensemble_info=validation.dict_to_JSlist(self.I.get_ensembles())
		else:
			ensemble_info=None
		Template_Dict['ensemble_info']=ensemble_info
		Template_Dict['sphere']=self.I.check_sphere()
		Template_Dict['num_ensembles']=self.I.check_ensembles()
		RB,flex,RB_nos,all_nos=self.I.get_RB_flex_dict()
		Template_Dict['Rigid_Body']=RB_nos
		Template_Dict['Flexible_Unit']=all_nos-RB_nos
		Template_Dict['RB_list']=validation.dict_to_JSlist_rows(RB,flex)
		Template_Dict['RB']=validation.utility.get_RB(validation.dict_to_JSlist_rows(RB,flex))
		Template_Dict['flex']=validation.utility.get_flex(validation.dict_to_JSlist_rows(RB,flex))
		Template_Dict['ID']=self.I.get_id()
		Template_Dict['ID_w']=self.I.get_id().split()
		Template_Dict['ID_T']=self.I.get_id()[0:6]+'_'+self.I.get_id()[6:]
		Template_Dict['ID_R']=(self.I.get_id()[0:6]+'_'+self.I.get_id()[6:]).split()
		Template_Dict['Molecule']=self.I.get_struc_title()
		Template_Dict['Title']=self.I.get_title()
		Template_Dict['Authors']=self.I.get_authors()
		Template_Dict['Entry_list']=validation.dict_to_JSlist(self.I.get_composition())
		Template_Dict['number_of_molecules']=self.I.get_number_of_models()
		Template_Dict['model_names']=self.I.get_model_names()
		Template_Dict['number_of_software']=self.I.get_software_length()
		Template_Dict['soft_list']=validation.dict_to_JSlist(self.I.get_software_comp())
		Template_Dict['number_of_datasets']=self.I.get_dataset_length()
		Template_Dict['Data']=[i.upper() for i in list(set(self.I.get_dataset_comp()['Dataset type']).difference({'Experimental model','Comparative model'}))]
		Template_Dict['Datasets_list']=validation.dict_to_JSlist(self.I.get_dataset_comp())
		Template_Dict['Protocols_number']=self.I.get_protocol_number()
		Template_Dict['Sampling_list']=validation.dict_to_JSlist(self.I.get_sampling())
		Template_Dict['num_chains']=int(len(self.I.get_composition()['Chain ID']))/int(len(list(Counter(self.I.get_composition()['Model ID']).keys())))
		return Template_Dict

	def run_model_quality(self,Template_Dict):
		print ("exo",self.I.check_sphere())
		if self.I.check_sphere()<1:
			#global clashscore; global rama; global sidechain;
			exv_data=None
			I_mp=get_molprobity_information.get_molprobity_information(self.mmcif_file)
			if I_mp.check_for_molprobity():
				dirname=os.path.normpath(os.getcwd() + os.sep + os.pardir)
				filename = os.path.abspath(os.path.join(dirname, 'working/Output/static/results/',str(Template_Dict['ID'])+'_temp_mp.txt'))
				if os.path.exists(filename):
					d_mp={}
					print ("Molprobity analysis file already exists...\n...assuming clashscores, Ramachandran and rotamer outliers have already been calculated")
					with open(filename,'rb') as fp:
						d_mp['molprobity']=pickle.load(fp)
					f_rota=os.path.abspath(os.path.join(dirname, 'working/Output/static/results/',str(Template_Dict['ID'])+'_temp_rota.txt'))
					with open(f_rota,'rb') as fp:
						d_mp['rota']=pickle.load(fp)
					f_rama=os.path.abspath(os.path.join(dirname, 'working/Output/static/results/',str(Template_Dict['ID'])+'_temp_rama.txt'))
					with open(f_rama,'rb') as fp:
						d_mp['rama']=pickle.load(fp)
					f_clash=os.path.abspath(os.path.join(dirname, 'working/Output/static/results/',str(Template_Dict['ID'])+'_temp_clash.txt'))
					with open(f_clash,'rb') as fp:
						d_mp['clash']=pickle.load(fp)
				else:
					print ("Molprobity analysis is being calculated...")
					manager = Manager()
					d_mp=manager.dict()
					validation.utility.runInParallel(I_mp.run_clashscore(d_mp),I_mp.run_ramalyze(d_mp),I_mp.run_rotalyze(d_mp),I_mp.run_molprobity(d_mp))
				a,b=I_mp.process_molprobity(d_mp['molprobity'])
				Template_Dict['bond']=len(a); Template_Dict['angle']=len(b)
				global clashscore;global rama;global sidechain
				clashscore,rama,sidechain=I_mp.get_data_for_quality_at_glance(d_mp['molprobity'])
				Template_Dict['molp_b']=validation.dict_to_JSlist(I_mp.molprobity_detailed_table_bonds(a))
				Template_Dict['molp_a']=validation.dict_to_JSlist(I_mp.molprobity_detailed_table_angles(b))
				Template_Dict['rotascore']=validation.dict_to_JSlist(I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
				Template_Dict['rotalist']=validation.dict_to_JSlist(I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota'])))
				Template_Dict['ramascore']=validation.dict_to_JSlist(I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
				Template_Dict['ramalist']=validation.dict_to_JSlist(I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama'])))
				clashscores,Template_Dict['tot']=I_mp.clash_summary_table(d_mp['clash'])
				Template_Dict['clashscore_list']=validation.dict_to_JSlist(clashscores)
				Template_Dict['clashlist']=I_mp.clash_detailed_table(d_mp['clash'])
				Template_Dict['assess_atomic_segments']='Clashscore: '+ str(clashscore) + ', Ramachandran outliers: '+ str(rama)+ ', Sidechain outliers: '+str(sidechain)
				Template_Dict['assess_excluded_volume']=['Not applicable']
			else:
				self.I.rewrite_mmcif()
				if I_mp.check_for_molprobity():
					print ("Molprobity analysis is being calculated...")
					manager = Manager()
					d_mp=manager.dict()
					runInParallel(I_mp.run_clashscore(d_mp),I_mp.run_ramalyze(d_mp),I_mp.run_rotalyze(d_mp),I_mp.run_molprobity(d_mp))
					a,b=I_mp.process_molprobity(d_mp['molprobity'])
					Template_Dict['bond']=len(a); Template_Dict['angle']=len(b)
					clashscore,rama,sidechain=I_mp.get_data_for_quality_at_glance(d_mp['molprobity'])
					Template_Dict['molp_b']=validation.dict_to_JSlist(I_mp.molprobity_detailed_table_bonds(a))
					Template_Dict['molp_a']=validation.dict_to_JSlist(I_mp.molprobity_detailed_table_angles(b))
					Template_Dict['rotascore']=validation.dict_to_JSlist(I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
					Template_Dict['rotalist']=validation.dict_to_JSlist(I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota'])))
					Template_Dict['ramascore']=validation.dict_to_JSlist(I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
					Template_Dict['ramalist']=validation.dict_to_JSlist(I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama'])))
					clashscores,Template_Dict['tot']=I_mp.clash_summary_table(d_mp['clash'])
					Template_Dict['clashscore_list']=validation.dict_to_JSlist(clashscores)
					Template_Dict['clashlist']=I_mp.clash_detailed_table(d_mp['clash'])
					Template_Dict['assess_atomic_segments']='Clashscore: '+ str(clashscore) + ', Ramachandran outliers: '+ str(rama)+ ', Sidechain outliers: '+str(sidechain)
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
				I_ev=get_excluded_volume.get_excluded_volume(self.mmcif_file)
				model_dict=I_ev.get_all_spheres()
				exv_data=I_ev.run_exc_vol_parallel(model_dict)
			
			Template_Dict['excluded_volume']=validation.dict_to_JSlist(exv_data)
			Template_Dict['assess_excluded_volume']=validation.utility.exv_readable_format(exv_data)
			clashscore=None
			rama=None
			sidechain=None
		return Template_Dict,clashscore,rama,sidechain,exv_data


	def run_sas_validation(self,Template_Dict):
		#global sas_data; global sas_fit;
		if self.I.check_for_sas():
			Template_Dict['sas']=["True"]
			I_sas=sas_validation.sas_validation(self.mmcif_file)
			Template_Dict['sasdb_code']=I_sas.get_SASBDB_code()
			Template_Dict['parameters_volume']=validation.dict_to_JSlist(I_sas.get_parameters_vol_many())
			Template_Dict['parameters_mw']=validation.dict_to_JSlist(I_sas.get_parameters_mw_many())
			Template_Dict['pddf_info']=validation.dict_to_JSlist(I_sas.get_pddf_info())
			print (Template_Dict['pddf_info'])
			Template_Dict['number_of_fits']=I_sas.get_total_fits()
			Template_Dict['chi_table']=validation.dict_to_JSlist(I_sas.get_chi_table())
			Template_Dict['rg_table']=validation.dict_to_JSlist(I_sas.get_rg_table_many())
			Template_Dict['sasdb_code_fits']=I_sas.get_sasdb_code_fits()
			I_sas_plt=validation.sas_validation_plots.sas_validation_plots(self.mmcif_file)
			I_sas.modify_intensity()
			I_sas_plt.plot_multiple()
			I_sas_plt.plot_pf()
			I_sas_plt.plot_Guinier()
			#if Template_Dict['number_of_fits']>0:
			I_sas_plt.plot_fits()
			#I_sas_plt.plot_residuals()
			I_sas.get_fit_image()
			sas_data=I_sas.get_rg_for_plot()
			sas_fit=I_sas.get_fits_for_plot()
		else:
			sas_data={}
			sas_fit={}
		return Template_Dict,sas_data,sas_fit

	def run_quality_glance(self,clashscore,rama,sidechain,exv_data,sas_data,sas_fit):
		I_plt=validation.get_plots.plots(self.mmcif_file)
		I_plt.plot_quality_at_glance(clashscore,rama,sidechain,exv_data,sas_data,sas_fit)

	def run_supplementary_table(self,args,Template_Dict):
		if (self.I.get_ensembles() is not None) and  (all_same(self.I.get_ensembles()['Clustering method'])):
			Template_Dict['clustering']=self.I.get_ensembles()['Clustering method'][0]
		elif self.I.get_ensembles() is not None:
			Template_Dict['clustering']=', '.join(self.I.get_ensembles()['Clustering method'])
		else:
			Template_Dict['clustering']='Not applicable'
		Template_Dict['location']=args.ls
		Template_Dict['complex_name']=self.I.get_struc_title()
		Template_Dict['PDB_ID']=self.I.get_id()
		Template_Dict['Subunits']=get_subunits(self.I.get_composition())
		Template_Dict['datasets']=get_datasets(self.I.get_dataset_details()) if self.I.get_dataset_details() is not None else 'Not provided or used'
		Template_Dict['physics']=physics
		Template_Dict['software']=get_software(self.I.get_software_comp())+args.ls
		Template_Dict['struc']=self.I.get_atomic_coverage()
		Template_Dict['method']=get_method_name(self.I.get_sampling())
		Template_Dict['method_type']=get_method_type(self.I.get_sampling())
		Template_Dict['method_details']=args.m
		Template_Dict['models']=', '.join(self.I.get_ensembles()['Number of models']) if self.I.get_ensembles() is not None else 'Not applicable' 
		Template_Dict['sampling_validation']=args.sv
		Template_Dict['feature']=self.I.get_ensembles()['Clustering feature'][0] if self.I.get_ensembles() is not None else 'Not applicable'
		Template_Dict['validation_input']=args.v1
		Template_Dict['cross_validation']=args.v2
		Template_Dict['model_precision']=', '.join([i+'&#8491' for i in self.I.get_ensembles()['Cluster precision']]) if self.I.get_ensembles() is not None else 'Model precision can not be calculated with one structure'
		Template_Dict['restraint_info']=get_restraints_info(self.I.get_restraints()) if self.I.get_restraints() is not None else 'Not provided or used'
		Template_Dict['Data_quality']=args.dv
		Template_Dict['clustering']=args.c
		Template_Dict['sampling_validation']=args.sv
		Template_Dict['resolution']=args.res
		return Template_Dict





