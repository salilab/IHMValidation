import pandas as pd
import glob
import sys,os,math
import numpy as np
import validation
import ihm
import ihm.reader
import re,pickle,requests,json
import multiprocessing as mp

class em_validation(validation.get_input_information):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(validation.get_input_information.get_id(self))
        self.nos=validation.get_input_information.get_number_of_models(self)
        self.dataset=validation.get_input_information.get_dataset_comp(self) 
    
