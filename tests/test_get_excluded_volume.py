import os,sys,glob
import unittest  
import pandas as pd
from io import StringIO, BytesIO
sys.path.insert(0, "../master/pyext/src/")
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
        check_xyz={1:[389.993,145.089,134.782,4.931],\
                    2:[406.895,142.176,135.653,3.318]}

        check_xyz_df = pd.DataFrame(data=check_xyz,index=['X','Y','Z','R'])

        model_dict=self.IO.get_all_spheres(filetemp=fh)
        list_of_sphere_list=list(model_dict.values())
        xyz_df=self.IO.get_xyzr(list_of_sphere_list[0])
        self.assertEqual(check_xyz_df.values.tolist(),xyz_df.values.tolist())

        add_chain={1:['A',1],2:['B',1]}
        add_chain_df = pd.DataFrame(data=add_chain,index=['Chain_ID','Model_ID'])

        fin=pd.concat([check_xyz_df,add_chain_df])
        xyz_complete_df=self.IO.get_xyzr_complete(model_ID=1,spheres=list_of_sphere_list[0])
        self.assertEqual(fin.values.tolist(),xyz_complete_df.values.tolist())

        viol_dict=self.IO.get_violation_dict(xyz_df)
        self.assertEqual({1: 0.0},self.IO.get_violation_dict(xyz_df))

        perc_satisfied=self.IO.get_violation_percentage(models_spheres_df=xyz_df,viols=viol_dict)
        self.assertEqual(100.0,perc_satisfied)


    def test_get_violation_others(self):
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
2 1 7 7 B 389.895 142.176 135.653 3.318 1.34 1
""")

        model_dict=self.IO.get_all_spheres(filetemp=fh)
        list_of_sphere_list=list(model_dict.values())
        xyz_df=self.IO.get_xyzr(list_of_sphere_list[0])

        viol_dict=self.IO.get_violation_dict(xyz_df)
        self.assertEqual({1: 1.0},viol_dict)

        perc_satisfied=self.IO.get_violation_percentage(models_spheres_df=xyz_df,viols=viol_dict)
        self.assertEqual(0.0,perc_satisfied)


    def test_violatio_multiple_models(self):
        fh = StringIO("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
2 . 1 1 1

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
2 1 7 7 B 389.895 142.176 135.653 3.318 1.34 1
3 1 1 6 A 489.993 145.089 134.782 4.931 0 2
4 1 7 7 B 589.895 142.176 135.653 3.318 1.34 2
""")

        model_dict=self.IO.get_all_spheres(filetemp=fh)
        output={'Models': [1, 2], 'Excluded Volume Satisfaction': [0.0, 100.0]}
        self.assertEqual(output,(self.IO.get_exc_vol_for_models_normalized(model_dict))

if __name__ == '__main__':
    unittest.main(warnings='ignore')

