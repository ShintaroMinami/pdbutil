pdbutil.py
===

A simple module for handling protein backbone coordinates.

## Requirement
* python3.x
* numpy
* biopython

## Usage
#### PDB file Read & Write
``` python
## Import functions
from pdbutil import read_pdb, write_pdb

## Read PDB file
data_dict = read_pdb('pdb_file_path.pdb')

data_dict: {
    'xyz_ca':    np.array [L, 3],    # C-alpha coordinates
    'xyz_bb':    np.array [L, 4, 3], # Backbone coordinates
    'xyz_aa':    np.array [L,14, 3], # All-atom coordinates
    'mask_aa':   np.array [L,14],    # Boolian array to address existing (=True) atoms
    'chain':     np.array [L,],      # Chain ID
    'resnum':    np.array [L,],      # Residue number
    'res1':      np.array [L,],      # One letter AA type
    'res3':      np.array [L,],      # Three letter AA type
    'occupancy': np.array [L,],      # Occupancy
    'bfactor':   np.array [L,],      # B-factor
    'insertion': np.array [L,],      # Insertion code
    'pdbstring': str,                # PDB format string
}

## Get PDB format string
pdb_string = write_pdb(**data_dict)   # Full writing
pdb_string = write_pdb(xyz_ca=xyz_ca) # Only xyz is required

# Write as PDB file
write_pdb(**data_dict, pdb_file='output.pdb')
```

#### Superpose & RMSD
```python
from pdbutil import superpose, calc_rmsd

## Superpose "Target" xyz (one or more) onto a "Reference" xyz
# xyz_reference.shape -> (L,3), (1,L,3), (L,4,3) or (1,L,4,3)
# xyz_targets.shape -> (L,3), (B,L,3), (L,4,3) or (B,L,4,3)

xyz_sup = superpose(xyz_reference, xyz_targets)

## Calculate C-alpha RMSD for all possible pairs
# xyz_ca1.shape -> (L,3), (1,L,3) or (B1,L,3)
# xyz_ca2.shape -> (L,3), (1,L,3) or (B2,L,3)

# Case of (L, 3) x (L, 3) -> float
rmsd = calc_rmsd(xyz_ca1, xyz_ca2)

# Case of (B1, L, 3) x (B2, L, 3) -> (B1, B2)
rmsd_matrix = calc_rmsd(xyz_ca1, xyz_ca2)
```

#### FASTA file Read & Write
``` python
## Import functions
from pdbutil import read_fasta

## Read FASTA file -> fasta object
fasta = read_fasta('fasta_file_path.fasta')

## Get deflines/sequences
deflines = fasta.deflines   # list[str], e.g. ['protein1', 'protein2', ...]
sequences = fasta.sequences # list[str], e.g. ['SEQUENCE1', 'SEQUENCE2', ...]

## Use the object as an Iterator
for data in fasta:
    print(data.defline)
    print(data.sequence)

## To get FASTA format string
fasta_string = str(fasta)

## To print as FASTA format on STDOUT
print(fasta)

```

## Author
* Shintaro Minami(https://github.com/ShintaroMinami)
* shintaro.minami@gmail.com

## License
MIT(https://choosealicense.com/licenses/mit/)

