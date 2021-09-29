import os
import sys
import unittest
import warnings
import tempfile
import ihm

path = os.path.abspath("master/pyext/src/")
sys.path.insert(0, path)
print (path)
from validation import GetInputInformation


def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            test_func(self, *args, **kwargs)
    return do_test


class Testing(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(Testing, self).__init__(*args, **kwargs)
        # self.mmcif_test_file='test.cif'
        # self.IO=GetInputInformation(self.mmcif_test_file)

    def test_basic_info(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""

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

                    """)
            self.IO = GetInputInformation(tmpfilepath)
            self.assertEqual('PDBDEVtest', self.IO.get_id())
            self.assertEqual('PDBDEVtest', self.IO.get_id_from_entry())
            self.assertEqual(
                'Title not available/Citation not provided', self.IO.get_title())
            self.assertEqual('Citation not present in file',
                             self.IO.get_authors())
            self.assertEqual(['1'], self.IO.get_assembly_ID_of_models())
            self.assertEqual(['1'], self.IO.get_representation_ID_of_models())

    def test_get_assembly_ID(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1 1 1
2 . 2 1 2
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(['1', '2'], temp.get_assembly_ID_of_models())

    def test_get_representation_ID(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 . 1A 1 AB
2 . 2B 1 BC
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(['1A', '2B'], temp.get_assembly_ID_of_models())
            self.assertEqual(
                ['AB', 'BC'], temp.get_representation_ID_of_models())

    def test_get_model_names(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1A 1 AB
2 Ialsolikeapple 2B 1 BC
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(['Ilikebanana', 'Ialsolikeapple'],
                             temp.get_model_names())

    def test_get_model_assem_dict(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1 1 AB
2 Ialsolikeapple 1 1 BC
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual({1: 1, 2: 1}, temp.get_model_assem_dict())

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1 1 AB
2 Ialsolikeapple X 1 BC
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertRaises(ValueError, lambda: temp.get_model_assem_dict())

    def test_get_model_rep_dict(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1 1 4
2 Ialsolikeapple 1 1 5 
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual({1: 4, 2: 5}, temp.get_model_rep_dict())

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1 1 AB
2 Ialsolikeapple X 1 BC

                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertRaises(ValueError, lambda: temp.get_model_rep_dict())

    def test_get_number_of_models(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 Ilikebanana 1 1 4
2 Ialsolikeapple 1 1 5
                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(2, temp.get_number_of_models())

    def test_get_protocol_number(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 'Best scoring model' 1 1 1

                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(1, temp.get_protocol_number())

    def test_get_sampling(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
1 'The best scoring model I could find, given available time management skills' 1 1 1
loop_
_ihm_modeling_protocol_details.id
_ihm_modeling_protocol_details.protocol_id
_ihm_modeling_protocol_details.step_id
_ihm_modeling_protocol_details.multi_scale_flag
_ihm_modeling_protocol_details.multi_state_flag
_ihm_modeling_protocol_details.step_name
_ihm_modeling_protocol_details.step_method
_ihm_modeling_protocol_details.num_models_end
1 1 1 YES NO banana apple 22
                    """)
            output = {'Step number': ['1'], 'Protocol ID': [1],
                      'Method name': ['apple'], 'Method type': ['banana'],
                      'Number of computed models': [22], 'Multi state modeling': ['False'],
                      'Multi scale modeling': ['True']}
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(output, temp.get_sampling())

    def test_get_chains_entities_assemblies(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
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

                    """)
            output = [('A', 'kiwi', 'kiwi', '1', 0),
                      ('B', 'berry', 'berry', '2', 1)]
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual([2], temp.get_number_of_chains())
            self.assertEqual(1, temp.get_number_of_assemblies())
            self.assertEqual(2, temp.get_number_of_entities())
            self.assertEqual(output, temp.get_all_asym())
            self.assertEqual({'A': [], 'B': []}, temp.get_empty_chain_dict())
            self.assertEqual({'A': 'kiwi', 'B': 'berry'},
                             temp.get_chain_subunit_dict())

    def test_get_rep(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""

loop_
_ihm_entity_poly_segment.id
_ihm_entity_poly_segment.entity_id
_ihm_entity_poly_segment.seq_id_begin
_ihm_entity_poly_segment.seq_id_end
1 1 1 6
#
loop_
_ihm_model_representation_details.id
_ihm_model_representation_details.representation_id
_ihm_model_representation_details.entity_id
_ihm_model_representation_details.entity_description
_ihm_model_representation_details.entity_asym_id
_ihm_model_representation_details.entity_poly_segment_id
_ihm_model_representation_details.model_object_primitive
_ihm_model_representation_details.starting_model_id
_ihm_model_representation_details.model_mode
_ihm_model_representation_details.model_granularity
_ihm_model_representation_details.model_object_count
1 1 1 apple A 1 sphere . rigid by-residue 1 


                    """)
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(
                ({'A': [['1-6', 'None']]}, {'A': []}, 1, 1), temp.get_RB_flex_dict())
            self.assertEqual('100%', temp.get_atomic_coverage())

    def test_get_software_info(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""

loop_
_software.pdbx_ordinal
_software.name
_software.classification
_software.description
_software.version
_software.type
_software.location
1 'haddock' 'docking'
'docking random things' developer program
https:chocolatemuffin.org

                    """)
            output = {'ID': ['1'], 
                    'Software name': ['<a href="https://pubmed.ncbi.nlm.nih.gov/12580598/">haddock</a>'], 
                    'Software version': ['developer'], 'Software classification': ['docking'], 
                    'Software location': ['<a href="https:chocolatemuffin.org">https:chocolatemuffin.org</a>']}
            temp = GetInputInformation(tmpfilepath)
            # print (temp.get_software_comp())
            self.assertEqual(output, temp.get_software_comp())

    def test_check_ensemble(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""

loop_
_ihm_ensemble_info.ensemble_id
_ihm_ensemble_info.ensemble_name
_ihm_ensemble_info.post_process_id
_ihm_ensemble_info.model_group_id
_ihm_ensemble_info.ensemble_clustering_method
_ihm_ensemble_info.ensemble_clustering_feature
_ihm_ensemble_info.num_ensemble_models
_ihm_ensemble_info.num_ensemble_models_deposited
_ihm_ensemble_info.ensemble_precision_value
_ihm_ensemble_info.ensemble_file_id
1 'Cluster 0' 1 1 . RMSD 1257 1 666 9

                    """)
            output = {'Ensemble number': ['1'],
                      'Ensemble name': ['Cluster 0'],
                      'Model ID': ['1'], 'Number of models': ['1257'],
                      'Clustering method': ['None'], 'Clustering feature': ['RMSD'],
                      'Cluster precision': ['666.0']}
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(output, temp.get_ensembles())

    def test_check_restraints(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""

loop_
_ihm_3dem_restraint.id
_ihm_3dem_restraint.dataset_list_id
_ihm_3dem_restraint.fitting_method
_ihm_3dem_restraint.fitting_method_citation_id
_ihm_3dem_restraint.struct_assembly_id
_ihm_3dem_restraint.number_of_gaussians
_ihm_3dem_restraint.model_id
_ihm_3dem_restraint.cross_correlation_coefficient
1 26 'Gaussian mixture models' 9 2 400 1 .
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
1 1 Nup120 3 17 LYS Nup120 3 412 LYS DSS 25
loop_
_ihm_2dem_class_average_restraint.id
_ihm_2dem_class_average_restraint.dataset_list_id
_ihm_2dem_class_average_restraint.number_raw_micrographs
_ihm_2dem_class_average_restraint.pixel_size_width
_ihm_2dem_class_average_restraint.pixel_size_height
_ihm_2dem_class_average_restraint.image_resolution
_ihm_2dem_class_average_restraint.image_segment_flag
_ihm_2dem_class_average_restraint.number_of_projections
_ihm_2dem_class_average_restraint.struct_assembly_id
_ihm_2dem_class_average_restraint.details
1 10 . 2.030 2.030 50 NO 1000 1 .
loop_
_ihm_sas_restraint.id
_ihm_sas_restraint.dataset_list_id
_ihm_sas_restraint.model_id
_ihm_sas_restraint.struct_assembly_id
_ihm_sas_restraint.profile_segment_flag
_ihm_sas_restraint.fitting_atom_type
_ihm_sas_restraint.fitting_method
_ihm_sas_restraint.fitting_state
_ihm_sas_restraint.radius_of_gyration
_ihm_sas_restraint.chi_value
_ihm_sas_restraint.details
1 18 1 2 NO 'Heavy atoms' FoXS Single 17.800 1.130 .

                    """)
            output = {'ID': [1, 2, 3, 4],
                      'Dataset ID': ['26', '25', '10', '18'],
                      'Restraint type': ['EM3DRestraint', 'CrossLinkRestraint', 'EM2DRestraint', 'SASRestraint'],
                      'Restraint info': ['Gaussian mixture models, 400', 'DSS, 1 cross-links',
                                         'Number of micrographs: None, Image resolution: 50.0',
                                         'Assembly name: None Fitting method: FoXS Multi-state: False']}
            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(output, temp.get_restraints())

    def test_check_dataset(self):
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
loop_
_ihm_dataset_list.id
_ihm_dataset_list.data_type
_ihm_dataset_list.database_hosted
_ihm_dataset_list.details
1 'Experimental model' YES .
4 'Comparative model' NO .
5 'CX-MS data' NO .
6 'EM raw micrographs' NO .
7 '2DEM class average' NO .
loop_
_ihm_3dem_restraint.id
_ihm_3dem_restraint.dataset_list_id
_ihm_3dem_restraint.fitting_method
_ihm_3dem_restraint.fitting_method_citation_id
_ihm_3dem_restraint.struct_assembly_id
_ihm_3dem_restraint.number_of_gaussians
_ihm_3dem_restraint.model_id
_ihm_3dem_restraint.cross_correlation_coefficient
1 6 'Gaussian mixture models' 9 2 400 1 .
loop_
_ihm_2dem_class_average_restraint.id
_ihm_2dem_class_average_restraint.dataset_list_id
_ihm_2dem_class_average_restraint.number_raw_micrographs
_ihm_2dem_class_average_restraint.pixel_size_width
_ihm_2dem_class_average_restraint.pixel_size_height
_ihm_2dem_class_average_restraint.image_resolution
_ihm_2dem_class_average_restraint.image_segment_flag
_ihm_2dem_class_average_restraint.number_of_projections
_ihm_2dem_class_average_restraint.struct_assembly_id
_ihm_2dem_class_average_restraint.details
1 7 . 2.030 2.030 50 NO 1000 1 .
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
1 1 Nup120 3 17 LYS Nup120 3 412 LYS DSS 5
loop_
_ihm_dataset_related_db_reference.id
_ihm_dataset_related_db_reference.dataset_list_id
_ihm_dataset_related_db_reference.db_name
_ihm_dataset_related_db_reference.accession_code
_ihm_dataset_related_db_reference.version
_ihm_dataset_related_db_reference.details
1 1 PDB 3JRO . .

                    """)
            output = {'ID': ['1', '4', '5', '6', '7'],
                      'Dataset type': ['Experimental model', 'Comparative model',
                                       'CX-MS data', 'EM raw micrographs', '2DEM class average'],
                      'Database name': ['PDB', '', '', '', ''],
                      'Details': ['PDB ID: 3JRO', 'template PDB ID: Not available',
                                  'Linker name and number of cross-links: DSS, 1 cross-links',
                                  'EMDB ID: Not available', 'EMDB ID: Not available']}

            temp = GetInputInformation(tmpfilepath)
            self.assertEqual(output, temp.get_dataset_details())
            self.assertEqual(False, temp.check_for_sas(
                temp.get_dataset_details()))
            self.assertEqual(
                True, (temp.check_for_cx(temp.get_dataset_details())))
            self.assertEqual(
                True, (temp.check_for_em(temp.get_dataset_details())))

    def test_mmcif_get_lists(self):
        """Test AtomSiteHandler"""
        with tempfile.TemporaryDirectory() as tempdir:
            tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

            with open(tmpfilepath, 'w') as tmpfile:
                tmpfile.write("""
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
            temp = GetInputInformation(tmpfilepath)
            a, b, c, d = temp.mmcif_get_lists()
            self.assertEqual([[], ['loop_'], ['_atom_site.group_PDB']], a)
            self.assertEqual([], d)
            self.assertEqual('_atom_site.B_iso_or_equiv', b[16])
            self.assertEqual('_atom_site.occupancy', b[13])
            self.assertEqual('1', c[19][-3])


if __name__ == '__main__':
    unittest.main(warnings='ignore')


