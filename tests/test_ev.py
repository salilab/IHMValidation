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
from validation.excludedvolume import get_excluded_volume
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
        self.IO=get_excluded_volume(self.mmcif_test_file)

    def test_get_all_spheres(self):
        """Test AtomSiteHandler"""
        fh = StringIO("""
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
1 1 1 6 A 389.993 145.089 134.782 4.931 0 1
2 1 7 7 B 406.895 142.176 135.653 3.318 1.34 1
""")
        self.assertEqual(1,len(list(self.IO.get_all_spheres(filetemp=fh).keys())))

    def test_get_XYZ(self):
        """Test AtomSiteHandler"""
        fh = StringIO("""
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
1 1 1 6 A 389.993 145.089 134.782 4.931 0 1
2 1 7 7 B 406.895 142.176 135.653 3.318 1.34 1
""")

        model_dict=self.IO.get_all_spheres(filetemp=fh)
        list_of_sphere_list=list(model_dict.values())
        xyz_df=self.IO.get_xyzr(list_of_sphere_list[0])
        self.assertEqual(406.895,xyz_df.iloc[0,1])


    def test_get_violation_dict(self):
        """Test AtomSiteHandler"""
        fh = StringIO("""
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
1 1 1 6 A 389.993 145.089 134.782 4.931 0 1
2 1 7 7 B 406.895 142.176 135.653 3.318 1.34 1
""")

        model_dict=self.IO.get_all_spheres(filetemp=fh)
        list_of_sphere_list=list(model_dict.values())
        xyz_df=self.IO.get_xyzr(list_of_sphere_list[0])
        self.assertEqual(0.0,list(self.IO.get_violation_dict(xyz_df).values())[0])

if __name__ == '__main__':
    unittest.main(warnings='ignore')

