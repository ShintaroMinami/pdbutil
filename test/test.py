#! /usr/bin/env python

import sys
sys.path.append('../')
from pdbutil import ProteinBackbone as pdb

struct = pdb(file='6dg5A.pdb')
struct.addHA(force=True)
struct.addH(force=True)
struct.addO(force=True)
struct.addH(force=True)
struct.calc_dihedral()
struct.calc_distmat()
struct.check_chainbreak()
struct.get_nearestN(5)

struct.printpdb(open('output.pdb', 'w'))