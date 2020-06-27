import os,sys
import unittest  
import ihm
import ihm.reader
import pandas as pd
import glob
import os
import numpy as np
import re
from collections import defaultdict
from io import StringIO, BytesIO
sys.path.append('../')
from validation import get_input_information,utility
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
        self.mmcif_test_file='test.cif'
        self.IO=get_input_information(self.mmcif_test_file)

    def test_get_databases(self):
        a=5  
        b=len(self.IO.get_databases())
        self.assertEqual(a,b)

    def test_basic_info(self):
        a='PDBDEVXX'
        self.assertEqual(a,self.IO.get_id())
        self.assertEqual(a,self.IO.get_id_from_entry())
        self.assertEqual('Characterizing something',self.IO.get_title())
        self.assertEqual('AA B;CC D',self.IO.get_authors())
        self.assertEqual('Structure of something cool',self.IO.get_struc_title())

    def test_model_details(self):
        """
        most of the functions below are based on ihm functions, 
        check https://github.com/ihmwg/python-ihm/blob/master/test for
        detailed testing
        """
        a=0
        self.assertEqual(a,self.IO.check_sphere())
        self.assertEqual(['1','1','1'],self.IO.get_assembly_ID_of_models())
        self.assertEqual(['1','1','1'],self.IO.get_representation_ID_of_models())
        self.assertEqual(['DD/Model 1', 'EE/Model 2', 'FF/Model 3'],self.IO.get_model_names())
        self.assertEqual({1: 1, 2: 1, 3: 1},self.IO.get_model_assem_dict())
        self.assertEqual(3,self.IO.get_number_of_models())
        self.assertEqual({1: 1, 2: 1, 3: 1},self.IO.get_model_rep_dict())
        self.assertEqual({'Model ID': [1, 1, 2, 2, 3, 3], 'Subunit number': [1, 2, 1, 2, 1, 2], 
                'Subunit ID': ['1', '1', '1', '1', '1', '1'], 
                'Subunit name': ['Ubiquitin', 'Ubiquitin', 'Ubiquitin', 'Ubiquitin', 'Ubiquitin', 'Ubiquitin'], 
                'Chain ID': ['A', 'B', 'A', 'B', 'A', 'B'], 
                'Total residues': [76, 76, 76, 76, 76, 76]},
                self.IO.get_composition())
        self.assertEqual(1,self.IO.get_protocol_number())
        self.assertEqual({'Step number': ['1'], 'Protocol ID': [1], 'Method name': [None], 
            'Method type': [None], 'Number of computed models': [None], 
            'Multi state modeling': ['True'], 'Multi scale modeling': ['None']},
            self.IO.get_sampling())
        self.assertEqual({'A': [], 'B': []},self.IO.get_RB_flex_dict()[0])
        self.assertEqual({'A': [['1-76']], 'B': [['1-76']]},self.IO.get_RB_flex_dict()[1])
        self.assertEqual(0,self.IO.get_RB_flex_dict()[2])
        self.assertEqual(2,self.IO.get_RB_flex_dict()[3])
        self.assertEqual(2,self.IO.get_number_of_chains()[0])
        self.assertEqual(('A', 'Ubiquitin', 'Ubiquitin', '1', 0),
                        self.IO.get_all_asym()[0])
        self.assertEqual({'A': 76, 'B': 76},
                        self.IO.get_residues_subunit_dict())
        self.assertEqual('Linker name and number of cross-links: EGS, 1 cross-links',
                            self.IO.get_dataset_xl_info('4'))
        self.assertEqual('SAS data/SASDCG7',self.IO.get_dataset_dict()['1'])
        self.assertEqual(0,self.IO.check_ensembles())
        self.assertEqual(['EGS, 1 cross-links', 'BS3, 1 cross-links', 'BS2G, 1 cross-links', 'DST, 1 cross-links'],
                        self.IO.get_restraints()['Restraint info'])
        self.assertEqual(['SAS data', 'Experimental model', 'Experimental model', 'CX-MS data', 'Single molecule FRET data'],
                        self.IO.get_dataset_comp()['Dataset type'])
        self.assertEqual(True,self.IO.check_for_sas(self.IO.get_dataset_comp()))
        self.assertEqual(True,self.IO.check_for_cx(self.IO.get_dataset_comp()))
        self.assertEqual(False,self.IO.check_for_em(self.IO.get_dataset_comp()))
        self.assertEqual('100%',self.IO.get_atomic_coverage())

    def test_mmcif_get_lists(self):
        """Test AtomSiteHandler"""
        fh = StringIO("""
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
ATOM 1 N N . SER 1 A 54.401 -49.984 -35.287 . 1 A 1 1 1
HETATM 2 C CA . SER . B 54.452 -48.492 -35.210 0.200 1 A 42.0 1 1
""")
        a,b,c,d=self.IO.mmcif_get_lists(filetemp=fh)
        self.assertEqual([[], ['loop_'], ['_atom_site.group_PDB']],a)
        self.assertEqual([],d)
        self.assertEqual('_atom_site.B_iso_or_equiv',b[16])
        self.assertEqual('_atom_site.occupancy',b[13])
        self.assertEqual('1',c[19][-3])


if __name__ == '__main__':
    unittest.main(warnings='ignore')

