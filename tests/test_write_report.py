import os
import sys
import unittest  
from collections import defaultdict
import warnings
import tempfile

os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'master', 'pyext', 'src'))
from validation import GetInputInformation, utility
from validation.Report import WriteReport

def ignore_warnings(test_func):
	def do_test(self, *args, **kwargs):
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", ResourceWarning)
			test_func(self, *args, **kwargs)
	return do_test

class Testing(unittest.TestCase):
	def __init__(self, *args, **kwargs):
		super(Testing, self).__init__(*args, **kwargs)
		self.template_dict=defaultdict()
		#self.mmcif_test_file='test.cif'
		#self.IO=cx_validation(self.mmcif_test_file)

	def test_run_entry_composition(self):
		with tempfile.TemporaryDirectory() as tempdir:
			tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

			with open(tmpfilepath, 'w') as tmpfile:
				tmpfile.write("""
data_PDBDEV_test
_entry.id PDBDEV_test
_struct.title 'Structure irrelevant'
loop_
_citation.id
_citation.title
1 "Structure of something cool"
loop_
_citation_author.citation_id
_citation_author.name
_citation_author.ordinal
1 'AA B' 1
1 'CC D' 2
loop_
_ihm_entity_poly_segment.id
_ihm_entity_poly_segment.entity_id
_ihm_entity_poly_segment.seq_id_begin
_ihm_entity_poly_segment.seq_id_end
1 1 1 6
2 1 7 20
#
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 kiwi
B 2 berry
#
loop_
_ihm_struct_assembly_details.id
_ihm_struct_assembly_details.assembly_id
_ihm_struct_assembly_details.parent_assembly_id
_ihm_struct_assembly_details.entity_description
_ihm_struct_assembly_details.entity_id
_ihm_struct_assembly_details.asym_id
_ihm_struct_assembly_details.entity_poly_segment_id
1 1 1 kiwi 1 A 1
2 1 1 berry 2 B 2
loop_
_entity.id
_entity.type
_entity.src_method
_entity.pdbx_description
_entity.formula_weight
_entity.pdbx_number_of_molecules
_entity.details
1 polymer man kiwi 96780.998 1 .
2 polymer man berry 98358.757 1 .
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1 1 AB
2 Ialsolikeapple 1 1 BC
					""")
			temp=WriteReport(tmpfilepath)
			d=temp.run_entry_composition(self.template_dict)
			self.assertEqual(None,d['ensemble_info'])
			self.assertEqual(0,d['sphere'])
			self.assertEqual(0,d['num_ensembles'])
			self.assertEqual(0,d['Rigid_Body'])
			self.assertEqual(0,d['Flexible_Unit'])
			self.assertEqual('PDBDEVtest',d['ID'])
			self.assertEqual(['PDBDEVtest'],d['ID_w'])
			self.assertEqual('PDBDEV_test',d['ID_T'])
			self.assertEqual(['PDBDEV_test'],d['ID_R'])
			self.assertEqual('Structure irrelevant',d['Molecule'])
			self.assertEqual('Structure of something cool',d['Title'])
			self.assertEqual('AA B; CC D',d['Authors'])
			self.assertEqual(2,d['number_of_molecules'])
			self.assertEqual(0,d['number_of_software'])
			self.assertEqual(['Ilikebanana', 'Ialsolikeapple'],d['model_names'])
			self.assertEqual([],d['soft_list'])
			self.assertEqual(0,d['number_of_datasets'])
			self.assertEqual(1,d['Protocols_number'])
			self.assertEqual(2.0,d['num_chains'])


if __name__ == '__main__':
    unittest.main(warnings='ignore')




