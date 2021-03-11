import os,sys
import unittest  
from io import StringIO, BytesIO
sys.path.insert(0, "../master/pyext/src/")
from validation import get_input_information,utility
import warnings,tempfile

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
