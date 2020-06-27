import os,sys
import unittest

sys.path.append("../../")
from FlaskApp import app
  
 
class BasicTests(unittest.TestCase):
 
    ############################
    #### setup and teardown ####
    ############################
 
    # executed prior to each test
    def setUp(self):
        app.config['TESTING'] = True
        app.config['DEBUG'] = False
        self.app = app.test_client()

    # executed after each test
    def tearDown(self):
        pass
 
 
###############
#### tests ####
###############
 
    def test_main_page(self):
        response = self.app.get('/', follow_redirects=True)
        self.assertEqual(response.status_code, 200)
 
    def test_main_page(self):
        response = self.app.get('/home', follow_redirects=True)
        self.assertEqual(response.status_code, 200)

 
if __name__ == "__main__":
    unittest.main()
