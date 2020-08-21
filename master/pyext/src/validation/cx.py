###################################
# Script : 
# 1) Contains class for XL-MS validation
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

from validation import get_input_information
import operator
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
from itertools import product

class cx_validation(get_input_information):
	def __init__(self,mmcif_file):
		super().__init__(mmcif_file)
		self.ID=str(get_input_information.get_id(self))
		self.nos=get_input_information.get_number_of_models(self)
		self.dataset=get_input_information.get_dataset_comp(self) 
	
	def get_xl_data (self):
		'''
		get a dataframe with information of XLs
		'''
		restraints = self.system.restraints
		xl_restraints=[r for r in restraints if 'cross' in str(r.__class__.__name__).lower()]
		lst=[]
		for r1 in xl_restraints:
			#linker=r1.linker._id
			linker_name=r1.linker.auth_name
			for xln in r1.experimental_cross_links:
				for xl in xln:
					res1_entity=xl.residue1.entity._id
					res2_entity=xl.residue2.entity._id
					res1_seq=xl.residue1.seq_id
					res2_seq=xl.residue2.seq_id
					res1_res2=str(res1_seq)+'_'+str(res2_seq)
					try:
						chains_1=self.get_asym_for_entity()[res1_entity]
						chains_2=self.get_asym_for_entity()[res2_entity]
						combinations=list(product(chains_1,chains_2))
						for i,j in enumerate(combinations):
							res1_asym=j[0]
							res2_asym=j[1]
							res1_id=str(res1_asym)+'_'+str(res1_seq)
							res2_id=str(res2_asym)+'_'+str(res2_seq)

							lst.append([linker_name,res1_res2,
								res1_entity,res1_seq,res1_asym,res1_id,
								res2_entity,res2_seq,res2_asym,res2_id])

					except: 
						pass

		xl_df=pd.DataFrame(lst,columns=['Linker_Name','XL_ID',
										'Res1_Entity_ID','Res1_Seq_ID','Res1_Chain','Res1_ID',
										'Res2_Entity_ID','Res2_Seq_ID','Res2_Chain','Res2_ID'])
		return xl_df

	def get_asym_for_entity(self):
		'''
		convert entity ID to chain ID
		'''
		asyms=self.system.asym_units
		asym_entity=dict()
		asym_entity={a.entity._id:[] for a in asyms}
		for a in asyms:
			asym_entity[a.entity._id].append(a._id)
		#asym_entity={a.entity._id:a._id for a in asyms}
		#asym_entity_again={a._id:a.entity._id for a in asyms}
		#print (asym_entity)
		return asym_entity

	def get_sphere_model_dict(self):
		'''
		get list of all spheres for all models present in mmcif 
		'''
		Model_object=[b for i in self.system.state_groups for j in i for a in j for b in a]
		model_dict={i+1:j._spheres for i,j in enumerate(Model_object)}
		return model_dict

	def get_atom_model_dict(self):
		'''
		get list of all atoms for all models present in mmcif 
		'''
		Model_object=[b for i in self.system.state_groups for j in i for a in j for b in a]
		model_dict={i+1:j._atoms for i,j in enumerate(Model_object)}
		return model_dict

	def get_xyzrseq_spheres(self,spheres):
		'''
		get xyzr,number of residues and information on structured/unstructured
		'''
		model_spheres={i+1:[j.seq_id_range,j.asym_unit._id,j.x,j.y,j.z,j.radius] for i,j in enumerate(spheres)}
		model_spheres_df=pd.DataFrame(model_spheres, index=['Seq','Chain','X','Y','Z','R']).T
		model_spheres_df['Res']=model_spheres_df['Seq'].apply(lambda x:int(x[1])-int(x[0])+1)
		model_spheres_df['Structured']=model_spheres_df['Res'].apply(lambda x: 0 if x>1 else 1)
		model_spheres_df_struc=model_spheres_df[model_spheres_df['Structured']==1]
		model_spheres_df_unstruc=model_spheres_df[model_spheres_df['Structured']==0]
		model_spheres_df_struc['Seq']=model_spheres_df_struc['Seq'].apply(lambda x: x[0])
		model_spheres_df_struc['Res_ID']=model_spheres_df_struc["Chain"]+'_' +model_spheres_df_struc["Seq"].astype(str)
		return model_spheres_df_struc,model_spheres_df_unstruc

	def get_xyzrseq_atoms(self,atoms):
		'''
		get xyz of CA atoms
		'''
		model_atoms={i+1:[j.seq_id,j.asym_unit._id,j.atom_id,j.x,j.y,j.z,] for i,j in enumerate(atoms)}
		model_atoms_df=pd.DataFrame(model_atoms, index=['Seq','Chain','Atom','X','Y','Z']).T
		model_atoms_df=model_atoms_df[model_atoms_df['Atom']=='CA']
		model_atoms_df['Res_ID']=model_atoms_df['Chain']+'_'+model_atoms_df['Seq'].astype(str)
		return model_atoms_df

	def convert_df_unstruc(seld,df):
		'''
		convert unstructured df into df that looks like a structured df
		'''
		lst=[]
		for index,row in df.iterrows():
			for j in range(row['Seq'][0],row['Seq'][1]+1):
				Res_ID=row['Chain']+'_'+str(j)
				lst.append([j,row['Chain'],row['X'],row['Y'],row['Z'],
							row['R'],row['Res'],row['Structured'],Res_ID])
		convert_df=pd.DataFrame(lst,columns=['Seq','Chain','X','Y','Z','R','Res','Structured','Res_ID'])
		return convert_df

	def get_complete_df_hybrid(self,xl_df,df):
		'''
		get complete df 
		add labels for struc, unstruc and hybrid xlinks
		extract XYZ of only xl residues
		'''
		lst=[]
		for index, row in xl_df.iterrows():
			try:
				df_res1=df[df['Res_ID']==row['Res1_ID']].iloc[:1,:] #.values.tolist()		
				df_res2=df[df['Res_ID']==row['Res2_ID']].iloc[:1,:] #.values.tolist()	
				struc_1=df_res1['Structured'].values[0];struc_2=df_res2['Structured'].values[0]
				if struc_1==struc_2 and struc_1==1:
					struc_value=1
				elif struc_1==struc_2 and struc_1==0:
					struc_value=0
				else:
					struc_value=2
				row_info=[row['Linker_Name'],row['XL_ID'],row['Res1_ID'],struc_value,
						df_res1['X'].values[0],df_res1['Y'].values[0],df_res1['Z'].values[0],row['Res2_ID'],
						df_res2['X'].values[0],df_res2['Y'].values[0],df_res2['Z'].values[0]]
				lst.append(row_info)
			except:
				pass	

		xl_df_comp=pd.DataFrame(lst,columns=['Linker','XL_ID','Res1','Structured','Res1_X','Res1_Y','Res1_Z',
									'Res2','Res2_X','Res2_Y','Res2_Z'])
		return xl_df_comp

	def get_complete_df_atomic(self,xl_df,df):
		'''
		get complete df 
		all labels are struc
		extract XYZ of only xl residues
		'''
		lst=[]
		for index,row in xl_df.iterrows():
			try:
				df_res1=df[df['Res_ID']==row['Res1_ID']].iloc[:1,:]
				df_res2=df[df['Res_ID']==row['Res2_ID']].iloc[:1,:]
				row_info=[row['Linker_Name'],row['XL_ID'],row['Res1_ID'],1,
				df_res1['X'].values[0],df_res1['Y'].values[0],df_res1['Z'].values[0],row['Res2_ID'],
				df_res2['X'].values[0],df_res2['Y'].values[0],df_res2['Z'].values[0]]
				lst.append(row_info)	
			except:
				pass
		xl_df_comp=pd.DataFrame(lst,columns=['Linker','XL_ID','Res1','Structured','Res1_X','Res1_Y','Res1_Z',
									'Res2','Res2_X','Res2_Y','Res2_Z'])
		return xl_df_comp

	def get_distance(self,df):
		'''
		distance between 2 residues
		'''
		df['dist']=((df['Res1_X']-df['Res2_X'])**2 +
					(df['Res1_Y']-df['Res2_Y'])**2 +
					(df['Res1_Z']-df['Res2_Z'])**2)**0.5
		return df


	def process_ambiguity(self,df):
		'''
		pick the smallest distance/xl if there are multiple values for the xl
		'''
		xl_list=list(df['dist'].groupby(df['XL_ID']))
		xl_dict={i[0]:sorted(i[1].values)[0] for i in xl_list if len(i[1].values)>1}
		xl_no_ambiguity=[i[0] for i in xl_list if len(i[1].values)==1]
		df_1=df[df['XL_ID'].isin(xl_no_ambiguity)]
		for key,val in xl_dict.items():
			df_2=df[(df['XL_ID']==key) & (df['dist']==val)]
			df_1=pd.concat([df_1,df_2])
		return df_1

	def label_inter_intra(self,df):
		'''
		label inter and intra differently
		'''
		df['Chain_A']=df['Res1'].apply(lambda x:x.split('_')[0])
		df['Chain_B']=df['Res2'].apply(lambda x:x.split('_')[0])
		df['Intra']=df.apply(lambda x:1 if x['Chain_A']==x['Chain_B'] else 0, axis=1)
		return df


	def get_violation(self,linker,dist):
		'''
		define violations based on linkers,
		needs to be updated with community standards
		'''
		if linker=='DSS' and dist<=30:
			return 1
		elif linker=='EDC' and dist<=20:
			return 1
		elif linker=='EDC' and dist>20:
			return 0
		elif dist<=30:
			return 1
		else:
			return 0

	def get_df_for_models(self,xl_df):
		'''
		get df for models 
		'''
		model_df=dict()
		if self.check_sphere()>0:
			model_dict=self.get_sphere_model_dict()
			for i, j in model_dict.items():
				df_struc,df_unstruc=self.get_xyzrseq_spheres(j)
				comp_df=self.convert_df_unstruc(df_unstruc)
				final_df=pd.concat((comp_df,df_struc),ignore_index = True)
				df_for_xl=self.get_complete_df_hybrid(xl_df,final_df)
				df_dist=self.get_distance(df_for_xl)
				df_intra=self.label_inter_intra(df_dist)
				df_intra['satisfaction']=df_intra.apply(lambda x: self.get_violation(x['Linker'], x['dist']), axis=1)
				df_final=self.process_ambiguity(df_intra)
				model_df[i]=df_final
		else:
			model_dict=self.get_atom_model_dict()
			for i, j in model_dict.items():
				df=self.get_xyzrseq_atoms(j)
				df_for_xl=self.get_complete_df_atomic(xl_df,df)
				df_dist=self.get_distance(df_for_xl)
				df_intra=self.label_inter_intra(df_dist)
				df_intra['satisfaction']=df_intra.apply(lambda x: self.get_violation(x['Linker'], x['dist']), axis=1)
				df_final=self.process_ambiguity(df_intra)
				model_df[i]=df_final
		return model_df


