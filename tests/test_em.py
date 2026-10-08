import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', 'ihm_validation')))
from em import EMValidation


# Trimmed from https://www.ebi.ac.uk/emdb/api/analysis/41983?information=all
# as EMDB served it with VA 0.0.1.dev137 (August 2026)
DEV137 = {
    '41983': {
        'version': '0.0.1.dev137',
        'resolution': {'value': 2.57},
        'recommended_contour_level': {'recl': 0.028, 'sigma': 2.97},
        'volume_estimate': {'volume': [30225.5444175397, 30225.5331380362],
                            'level': [-0.0954404846, -0.0925947464],
                            'estvolume': 485.92},
        'central_slice': {'indices': {'x': 175, 'y': 175, 'z': 175}},
        'feature_assessment': {'feature_zones': {'feature_zone': 0.0637}},
        'validation_json': {'status': False, 'schema_version': '0.0.1.dev137'},
    }
}

# The same entry as EMDB serves it with VA 0.0.1.dev151 (October 2026):
# the dev137 analysis moved one level down, under the EMDB ID
DEV151 = {
    '41983': {
        'version': '0.0.1.dev151',
        'resolution': {'value': 2.57},
        'recommended_contour_level': {'recl': 0.028, 'sigma': 2.97},
        'feature_assessment': {},
        'validation_json': {'status': False, 'schema_version': '0.0.1.dev151'},
        'audit_validation': {'cif_revision_audit': {'decision': 'run'}},
        'map_map_cc': {'emd_41983.map_emd_41983_half_map_1.map': {
            'names': ['emd_41983.map', 'emd_41983_half_map_1.map'],
            'value': 0.85}},
        '41983': DEV137['41983'],
    }
}


class TestUnwrapEmdbAnalysis(unittest.TestCase):

    def test_flat_analysis(self):
        out = EMValidation.unwrap_emdb_analysis(DEV137, '41983')
        self.assertEqual(out['volume_estimate']['estvolume'], 485.92)
        self.assertEqual(out['central_slice']['indices']['x'], 175)
        self.assertEqual(out['version'], '0.0.1.dev137')

    def test_nested_analysis(self):
        out = EMValidation.unwrap_emdb_analysis(DEV151, '41983')
        self.assertEqual(out['volume_estimate']['estvolume'], 485.92)
        self.assertEqual(out['central_slice']['indices']['x'], 175)
        self.assertEqual(out['recommended_contour_level']['recl'], 0.028)
        self.assertNotIn('41983', out)

    def test_nested_analysis_prefers_inner_values(self):
        # The outer copies describe the dev151 wrapper, not the analysis:
        # its feature assessment is empty
        out = EMValidation.unwrap_emdb_analysis(DEV151, '41983')
        self.assertEqual(out['version'], '0.0.1.dev137')
        self.assertEqual(out['feature_assessment'],
                         {'feature_zones': {'feature_zone': 0.0637}})

    def test_nested_analysis_keeps_outer_only_keys(self):
        out = EMValidation.unwrap_emdb_analysis(DEV151, '41983')
        self.assertEqual(
            out['map_map_cc']['emd_41983.map_emd_41983_half_map_1.map']['value'],
            0.85)
        self.assertEqual(
            out['audit_validation']['cif_revision_audit']['decision'], 'run')


if __name__ == '__main__':
    unittest.main()
