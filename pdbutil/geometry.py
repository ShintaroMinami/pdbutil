import numpy as np

_atom2id = {'N':0, 'CA':1, 'C':2, 'O':3, 'CB':4}
_param = {
    'angle_N_CA_C': np.deg2rad(111.2), 'angle_CA_C_N':np.deg2rad(116.2), 'angle_C_N_CA': np.deg2rad(121.7),
    'angle_N_CA_CB': np.deg2rad(110.6), 'angle_CB_CA_C': np.deg2rad(110.6), 'angle_C_N_H': np.deg2rad(123.0),
    'angle_N_C_O': np.deg2rad(122.7), 'angle_N_CA_1HA': np.deg2rad(109.5), 'angle_N_CA_2HA': np.deg2rad(109.5),
    'dhdrl_C_N_CA_CB': np.deg2rad(-121.4), 'dhdrl_N_C_CA_CB': np.deg2rad(121.4), 'dhdrl_CA_C_N_H': np.deg2rad(0.0),
    'dhdrl_CA_N_C_O': np.deg2rad(0.0), 'dhdrl_C_N_CA_1HA': np.deg2rad(121.4), 'dhdrl_C_N_CA_2HA': np.deg2rad(-121.4),
    'length_CN': 1.33, 'length_NCA': 1.46, 'length_CAC': 1.52,
    'length_CC': 1.54, 'length_CO': 1.24, 'length_NH': 1.01,
    'length_CH': 1.09
    }


## add virtual O atoms ##
def addO(xyz, force=False):
    xyz_O = []
    for iaa in range(len(xyz)-1):
        co = _zmat2xyz(
            _param['length_CO'],
            _param['angle_N_C_O'],
            _param['dhdrl_CA_N_C_O'],
            xyz[iaa+1][_atom2id['CA']],
            xyz[iaa+1][_atom2id['N']],
            xyz[iaa][_atom2id['C']]
            )
        xyz_O.append(co)
    return np.stack(xyz_O)

## add virtual CB atoms ##
'''
def addCB(self, force=False):
    for iaa in range(len(self.coord)):
        if ((self.exists[iaa][self.atom2id['CB']] == True) and (force==False)): continue
        cb1 = _zmat2xyz(self.param['length_CC'],
                        self.param['angle_N_CA_CB'],
                        self.param['dhdrl_C_N_CA_CB'],
                        self.coord[iaa][self.atom2id['C']],
                        self.coord[iaa][self.atom2id['N']],
                        self.coord[iaa][self.atom2id['CA']])
        cb2 = _zmat2xyz(self.param['length_CC'],
                        self.param['angle_CB_CA_C'],
                        self.param['dhdrl_N_C_CA_CB'],
                        self.coord[iaa][self.atom2id['N']],
                        self.coord[iaa][self.atom2id['C']],
                        self.coord[iaa][self.atom2id['CA']])
        cb = (cb1 + cb2)/2.0
        self.coord[iaa][self.atom2id['CB']][0] = cb[0]
        self.coord[iaa][self.atom2id['CB']][1] = cb[1]
        self.coord[iaa][self.atom2id['CB']][2] = cb[2]
        self.exists[iaa][self.atom2id['CB']] = True
'''







#### Functions ####
def _zmat2xyz(bond, angle, dihedral, one, two , three):
    oldvec = np.ones(4, dtype=np.float)
    oldvec[0] = bond * np.sin(angle) * np.sin(dihedral)
    oldvec[1] = bond * np.sin(angle) * np.cos(dihedral)
    oldvec[2] = bond * np.cos(angle)
    mat = _viewat(three, two, one)
    newvec = np.dot(mat, oldvec)
    # return
    return newvec

def _viewat(p1, p2, p3):
    # vector #
    p12 = p2 - p1
    p13 = p3 - p1
    # normalize #
    z = p12 / np.linalg.norm(p12)
    # crossproduct #
    x = np.cross(p13, p12)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    y /= np.linalg.norm(y)
    # transpation matrix
    mat = np.zeros((4, 4), dtype=np.float)
    for i in range(3):
        mat[i][0] = x[i]
        mat[i][1] = y[i]
        mat[i][2] = z[i]
        mat[i][3] = p1[i]
    mat[3][3] = 1.0
    # return
    return mat

def _xyz2dihedral(p1, p2, p3, p4):
    # small val #
    eps = 0.0000001
    # bond vector
    v1 = p2 - p1
    v2 = p3 - p2
    v3 = p4 - p3
    # perpendicular vector #
    perp123 = np.cross(v1, v2)
    perp234 = np.cross(v2, v3)
    perp123 /= np.linalg.norm(perp123)
    perp234 /= np.linalg.norm(perp234)
    # scalar product #
    scp = np.dot(perp123, perp234)
    scp = scp - eps if scp > 1 else scp
    scp = scp + eps if scp < -1 else scp
    # absolute angle #
    angle = np.rad2deg( np.arccos(scp) )
    # return #
    return angle if np.dot(v1, perp234) > 0 else -angle
