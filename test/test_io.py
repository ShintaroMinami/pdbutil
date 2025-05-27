#! /usr/bin/env python

import sys
sys.path.append('../')

from pdbutil import read_pdb, write_pdb

data_dict = read_pdb('6dg5A.pdb')
print(write_pdb(**data_dict, model_type='aa'))