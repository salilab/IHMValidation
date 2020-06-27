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
from validation.molprobity import get_molprobity_information
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
        self.IO=get_molprobity_information(self.mmcif_test_file)

    def test_molprobity_true(self):
        fh = StringIO("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 'Best scoring model' 1 2 3
2 'Best scoring model' 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
2 "Cluster 2" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
2 2
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
        self.assertEqual(True,self.IO.check_for_molprobity(filetemp=fh))

    def test_molprobity_false(self):
        fh = StringIO("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 'Best scoring model' 1 2 3
2 'Best scoring model' 1 1 1
#
loop_
_ihm_model_group.id
_ihm_model_group.name
_ihm_model_group.details
1 "Cluster 1" .
2 "Cluster 2" .
#
loop_
_ihm_model_group_link.group_id
_ihm_model_group_link.model_id
1 1
2 2
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
""")
        self.assertEqual(False,self.IO.check_for_molprobity(filetemp=fh))

    

if __name__ == '__main__':
    unittest.main(warnings='ignore')

