import os,sys
import unittest  
import ihm
import ihm.reader
import pandas as pd
import glob
import os,shutil
import numpy as np
import re
from collections import defaultdict
from io import StringIO, BytesIO
sys.path.insert(0, "../master/pyext/src/")
from validation import get_input_information,utility
from validation.cx import cx_validation
import warnings

def ignore_warnings(test_func):
	def do_test(self, *args, **kwargs):
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", ResourceWarning)
			test_func(self, *args, **kwargs)
	return do_test

class Testing(unittest.TestCase):
	def __init__(self, *args, **kwargs):
		super(Testing, self).__init__(*args, **kwargs)
		#self.mmcif_test_file='test.cif'
		#self.IO=cx_validation(self.mmcif_test_file)

	def test_get_xl_data_1(self):
		"""Test AtomSiteHandler"""
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_cross_link_list.id
_ihm_cross_link_list.group_id
_ihm_cross_link_list.entity_description_1
_ihm_cross_link_list.entity_id_1
_ihm_cross_link_list.seq_id_1
_ihm_cross_link_list.comp_id_1
_ihm_cross_link_list.entity_description_2
_ihm_cross_link_list.entity_id_2
_ihm_cross_link_list.seq_id_2
_ihm_cross_link_list.comp_id_2
_ihm_cross_link_list.linker_type
_ihm_cross_link_list.dataset_list_id
1 1 foo 1 2 THR foo 1 3 CYS DSS 97
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')

		self.assertEqual(1,I.get_xl_data().shape[0])
		self.assertEqual(10,I.get_xl_data().shape[1])
		self.assertEqual('DSS',I.get_xl_data()['Linker_Name'].unique()[0])
		self.assertEqual('1',I.get_xl_data()['Res1_Entity_ID'].unique()[0])
		self.assertEqual(2,I.get_xl_data()['Res1_Seq_ID'].unique()[0])
		self.assertEqual('1',I.get_xl_data()['Res2_Entity_ID'].unique()[0])
		self.assertEqual(3,I.get_xl_data()['Res2_Seq_ID'].unique()[0])
		self.assertEqual('A',I.get_xl_data()['Res1_Chain'].unique()[0])
		self.assertEqual('A',I.get_xl_data()['Res2_Chain'].unique()[0])


	def test_get_xl_data_2(self):
		"""Test AtomSiteHandler"""
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_cross_link_list.id
_ihm_cross_link_list.group_id
_ihm_cross_link_list.entity_description_1
_ihm_cross_link_list.entity_id_1
_ihm_cross_link_list.seq_id_1
_ihm_cross_link_list.comp_id_1
_ihm_cross_link_list.entity_description_2
_ihm_cross_link_list.entity_id_2
_ihm_cross_link_list.seq_id_2
_ihm_cross_link_list.comp_id_2
_ihm_cross_link_list.linker_type
_ihm_cross_link_list.dataset_list_id
1 1 foo 1 2 THR foo 1 3 CYS DSS 97
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		self.assertEqual('DSS',I.get_xl_data()['Linker_Name'].unique()[0])
		self.assertEqual('1',I.get_xl_data()['Res1_Entity_ID'].unique()[0])
		self.assertEqual(2,I.get_xl_data()['Res1_Seq_ID'].unique()[0])
		self.assertEqual('1',I.get_xl_data()['Res2_Entity_ID'].unique()[0])
		self.assertEqual(3,I.get_xl_data()['Res2_Seq_ID'].unique()[0])
		self.assertEqual('A',I.get_xl_data()['Res1_Chain'].unique()[0])
		self.assertEqual('A',I.get_xl_data()['Res2_Chain'].unique()[0])



	def test_get_asym_for_entity(self):
		"""Test AtomSiteHandler"""
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		self.assertEqual('1',list(I.get_asym_for_entity().keys())[0]) 
		self.assertEqual(['A'],list(I.get_asym_for_entity().values())[0]) 

	def test_get_asym_for_entity_dimer(self):
		"""Test AtomSiteHandler"""
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo
B 1 foo
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		self.assertEqual('1',list(I.get_asym_for_entity().keys())[0]) 
		self.assertEqual(['A','B'],list(I.get_asym_for_entity().values())[0]) 


	def test_get_sphere_model_dict(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_ihm_sphere_obj_site.id
_ihm_sphere_obj_site.entity_id
_ihm_sphere_obj_site.seq_id_begin
_ihm_sphere_obj_site.seq_id_end
_ihm_sphere_obj_site.asym_id
_ihm_sphere_obj_site.Cartn_x
_ihm_sphere_obj_site.Cartn_y
_ihm_sphere_obj_site.Cartn_z
_ihm_sphere_obj_site.object_radius
_ihm_sphere_obj_site.rmsf
_ihm_sphere_obj_site.model_id
1 1 1 6 A 389.993 145.089 134.782 4.931 . 1
2 1 7 7 B 406.895 142.176 135.653 3.318 1.34 1

"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		self.assertEqual(2,len(I.get_sphere_model_dict()[1]) )

	def test_get_atom_model_dict(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.label_entity_id
_atom_site.auth_asym_id
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
_atom_site.ihm_model_id
ATOM 1 N N . SER 1 A 54.401 -49.984 -35.287 . 1 A . 1 1
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		self.assertEqual(1,len(I.get_atom_model_dict()[1]) )

	def test_get_xyzrseq_spheres(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_ihm_sphere_obj_site.id
_ihm_sphere_obj_site.entity_id
_ihm_sphere_obj_site.seq_id_begin
_ihm_sphere_obj_site.seq_id_end
_ihm_sphere_obj_site.asym_id
_ihm_sphere_obj_site.Cartn_x
_ihm_sphere_obj_site.Cartn_y
_ihm_sphere_obj_site.Cartn_z
_ihm_sphere_obj_site.object_radius
_ihm_sphere_obj_site.rmsf
_ihm_sphere_obj_site.model_id
1 1 1 6 A 389.993 145.089 134.782 4.931 . 1
1 1 7 7 A 389.993 145.089 134.782 4.931 . 1

"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		spheres=I.get_sphere_model_dict()[1]
		self.assertEqual((1,9),I.get_xyzrseq_spheres(spheres)[0].shape)
		self.assertEqual((1,8),I.get_xyzrseq_spheres(spheres)[1].shape)
		self.assertEqual('A',I.get_xyzrseq_spheres(spheres)[0]['Chain'].unique()[0])
		self.assertEqual('A',I.get_xyzrseq_spheres(spheres)[1]['Chain'].unique()[0])
		self.assertEqual(389.993,I.get_xyzrseq_spheres(spheres)[0]['X'].unique()[0])
		self.assertEqual(145.089,I.get_xyzrseq_spheres(spheres)[1]['Y'].unique()[0])

	def test_get_atom_model_dict(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.label_entity_id
_atom_site.auth_asym_id
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
_atom_site.ihm_model_id
ATOM 1 CA CA . SER 1 A 54.401 -49.984 -35.287 . 1 A . 1 1
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		atom=I.get_atom_model_dict()[1]
		self.assertEqual(1,I.get_xyzrseq_atoms(atom)['Seq'].unique()[0])
		self.assertEqual('A',I.get_xyzrseq_atoms(atom)['Chain'].unique()[0])
		self.assertEqual('CA',I.get_xyzrseq_atoms(atom)['Atom'].unique()[0])
		self.assertEqual(54.401,I.get_xyzrseq_atoms(atom)['X'].unique()[0])
		self.assertEqual('A_1',I.get_xyzrseq_atoms(atom)['Res_ID'].unique()[0])

	def test_convert_df_unstruc(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_ihm_sphere_obj_site.id
_ihm_sphere_obj_site.entity_id
_ihm_sphere_obj_site.seq_id_begin
_ihm_sphere_obj_site.seq_id_end
_ihm_sphere_obj_site.asym_id
_ihm_sphere_obj_site.Cartn_x
_ihm_sphere_obj_site.Cartn_y
_ihm_sphere_obj_site.Cartn_z
_ihm_sphere_obj_site.object_radius
_ihm_sphere_obj_site.rmsf
_ihm_sphere_obj_site.model_id
1 1 1 6 A 389.993 145.089 134.782 4.931 . 1
1 1 7 7 A 389.993 145.089 134.782 4.931 . 1

"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		spheres=I.get_sphere_model_dict()[1]
		self.assertEqual((6,9),I.convert_df_unstruc(I.get_xyzrseq_spheres(spheres)[1]).shape)

	def test_get_complete_df_hybrid(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_ihm_sphere_obj_site.id
_ihm_sphere_obj_site.entity_id
_ihm_sphere_obj_site.seq_id_begin
_ihm_sphere_obj_site.seq_id_end
_ihm_sphere_obj_site.asym_id
_ihm_sphere_obj_site.Cartn_x
_ihm_sphere_obj_site.Cartn_y
_ihm_sphere_obj_site.Cartn_z
_ihm_sphere_obj_site.object_radius
_ihm_sphere_obj_site.rmsf
_ihm_sphere_obj_site.model_id
1 1 1 6 A 389.993 145.089 134.782 4.931 . 1
1 1 7 7 A 389.993 145.089 134.782 4.931 . 1
loop_
_ihm_cross_link_list.id
_ihm_cross_link_list.group_id
_ihm_cross_link_list.entity_description_1
_ihm_cross_link_list.entity_id_1
_ihm_cross_link_list.seq_id_1
_ihm_cross_link_list.comp_id_1
_ihm_cross_link_list.entity_description_2
_ihm_cross_link_list.entity_id_2
_ihm_cross_link_list.seq_id_2
_ihm_cross_link_list.comp_id_2
_ihm_cross_link_list.linker_type
_ihm_cross_link_list.dataset_list_id
1 1 foo 1 2 THR foo 1 3 CYS DSS 97
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo
"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		spheres=I.get_sphere_model_dict()[1]
		df_s,df_u=I.get_xyzrseq_spheres(spheres)
		df=I.convert_df_unstruc(df_u)
		xl_df=I.get_xl_data()
		self.assertEqual((1,11),I.get_complete_df_hybrid(xl_df,df).shape)
		self.assertEqual(389.993,I.get_complete_df_hybrid(xl_df,df)['Res1_X'].unique()[0])
		self.assertEqual('2_3',I.get_complete_df_hybrid(xl_df,df)['XL_ID'].unique()[0])
		self.assertEqual('A_2',I.get_complete_df_hybrid(xl_df,df)['Res1'].unique()[0])

	def test_get_atom_model_dict(self):
		fh = """
data_PDBDEV_test
_entry.id PDBDEV_test
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.label_entity_id
_atom_site.auth_asym_id
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
_atom_site.ihm_model_id
ATOM 1 CA CA . SER 1 A 54.401 -49.984 -35.287 . 1 A . 1 1
ATOM 2 CA CA . SER 2 A 54.401 -49.984 -35.287 . 3 A . 1 1
ATOM 3 CA CA . SER 3 A 54.401 -49.984 -35.287 . 3 A . 1 1

loop_
_ihm_cross_link_list.id
_ihm_cross_link_list.group_id
_ihm_cross_link_list.entity_description_1
_ihm_cross_link_list.entity_id_1
_ihm_cross_link_list.seq_id_1
_ihm_cross_link_list.comp_id_1
_ihm_cross_link_list.entity_description_2
_ihm_cross_link_list.entity_id_2
_ihm_cross_link_list.seq_id_2
_ihm_cross_link_list.comp_id_2
_ihm_cross_link_list.linker_type
_ihm_cross_link_list.dataset_list_id
1 1 foo 1 2 THR foo 1 3 CYS DSS 97
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo

"""    
		with open ('test.cif', 'w') as fd:
			fd.write(fh)
		I=cx_validation('test.cif')
		atom=I.get_atom_model_dict()[1]
		df=I.get_xyzrseq_atoms(atom)
		xl_df=I.get_xl_data()
		self.assertEqual((1,11),I.get_complete_df_atomic(xl_df,df).shape)
		self.assertEqual(54.401,I.get_complete_df_atomic(xl_df,df)['Res1_X'].unique()[0])
		self.assertEqual('2_3',I.get_complete_df_atomic(xl_df,df)['XL_ID'].unique()[0])
		self.assertEqual('A_2',I.get_complete_df_atomic(xl_df,df)['Res1'].unique()[0])

	def test_get_distance(self):
		lst=pd.DataFrame([[4,5,6,1,2,3]],columns=['Res1_X','Res1_Y','Res1_Z',
													'Res2_X','Res2_Y','Res2_Z'])
		I=cx_validation('test.cif')
		self.assertEqual(5,int(I.get_distance(lst)['dist'].values[0]))

	def test_label_inter_intra_1(self):
		lst=pd.DataFrame([['C_4','C_5']],columns=['Res1','Res2'])
		I=cx_validation('test.cif')
		self.assertEqual(1,int(I.label_inter_intra(lst)['Intra'].values[0]))

	def test_label_inter_intra_2(self):
		lst=pd.DataFrame([['C_4','D_5']],columns=['Res1','Res2'])
		I=cx_validation('test.cif')
		self.assertEqual(0,int(I.label_inter_intra(lst)['Intra'].values[0]))

	def test_get_violation(self):
		I=cx_validation('test.cif')
		print (I.get_violation('DSS',31))
		self.assertEqual(0,I.get_violation('DSS',31))
		self.assertEqual(1,I.get_violation('DSS',29))
		self.assertEqual(1,I.get_violation('EDC',20))
		self.assertEqual(0,I.get_violation('EDC',29))

	def test_process_ambiguity(self):
		lst=pd.DataFrame([['C_4',23],['C_4',11],['C_5',22]],columns=['XL_ID','dist'])
		I=cx_validation('test.cif')
		self.assertEqual(22,I.process_ambiguity(lst)['dist'].values[0])
		self.assertEqual(11,I.process_ambiguity(lst)['dist'].values[1])
		self.assertEqual('C_5',I.process_ambiguity(lst)['XL_ID'].values[0])
		self.assertEqual('C_4',I.process_ambiguity(lst)['XL_ID'].values[1])


if __name__ == '__main__':
	unittest.main(warnings='ignore')

