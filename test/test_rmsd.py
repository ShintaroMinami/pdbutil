#! /usr/bin/env python

import sys
sys.path.append('../')

import numpy as np
from scipy.spatial.transform import Rotation
from pdbutil import calc_rmsd, read_pdb, write_pdb
from pdbutil.rmsd import superpose


data_dict = read_pdb('6dg5A.pdb')
xyz = data_dict['xyz_bb']
xyz_aa = data_dict['xyz_aa']
print(data_dict['xyz_aa'].shape)

xyz2 = []
for _ in range(20):
    rot_random = Rotation.random().as_matrix()
    xyz2.append(np.einsum('i j, l a i -> l a j', rot_random, xyz) + np.random.randn(1,1,3))
xyz2 = np.stack(xyz2, axis=0)

rmsd = calc_rmsd(xyz2[:,:,2], xyz2[:10,:,1])
print(rmsd.shape)

sup, rmsd = superpose(xyz, xyz2, return_rmsd=True)
data_dict['xyz_bb'] = sup[0]

print(rmsd)

