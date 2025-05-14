
def pdb_rmsd():
    """
    Command line interface for calculating RMSD.
    """
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Calculate RMSD between two PDB files.', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('pdb1', type=str, nargs='+', help='First groupd of PDB files.')
    parser.add_argument('--pdb2', type=str, nargs='+', help='Second group of PDB files.')
    parser.add_argument('--csv', action='store_true', help='Output in CSV format.')
    args = parser.parse_args()

    import numpy as np
    from pdbutil.pdbio import read_pdb    
    from pdbutil.rmsd import calc_rmsd
    
    xyz1, file1 = [read_pdb(f)['xyz_ca'] for f in args.pdb1], args.pdb1

    if not np.all(np.array([len(x) == len(xyz1[0]) for x in xyz1], dtype=bool)):
        raise ValueError("All input PDB files must have the same number of CA atoms.")

    (xyz2, file2) = ([read_pdb(f)['xyz_ca'] for f in args.pdb2], args.pdb2) if args.pdb2 is not None else (xyz1, file1)

    if not np.all(np.array([len(x) == len(xyz1[0]) for x in xyz2], dtype=bool)):
        raise ValueError("All input PDB files must have the same number of CA atoms.")

    rmsd_mat = calc_rmsd(np.stack(xyz1, axis=0), np.stack(xyz2, axis=0))
    file_mat = np.array(file1)[:,None] + ":" + np.array(file2)[None,:]
    if args.csv:
        print("RMSD,File1,File2")
    else:
        print(f" {'RMSD':>8}  {'File1':<20}  {'File2':<20}")
    for rmsd, file_pair in zip(rmsd_mat.flatten(), file_mat.flatten()):
        if args.csv:
            vals = [f"{rmsd:.5f}"] + file_pair.split(':')
            print(",".join(vals))
        else:
            print(f" {rmsd:8.5f}  {file_pair}")


def pdb_superpose():
    """
    Command line interface for superposing PDB structures onto a reference.
    """
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Superpose two PDB files.', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('pdbs', type=str, nargs='+', help='Target PDB files.')
    parser.add_argument('-r', '--reference', type=str, default=None, help='Reference PDB file.')
    parser.add_argument('-o', '--output_dir', type=str, default="", help='Output directory.')
    args = parser.parse_args()

    from pathlib import Path
    import numpy as np
    from pdbutil.pdbio import read_pdb, write_pdb
    from pdbutil.rmsd import superpose
    
    datadicts, names = [read_pdb(f) for f in args.pdbs], [str(Path(f).stem) for f in args.pdbs]
    
    xyz_ref = datadicts[0]['xyz_bb'] if args.reference is None else read_pdb(args.reference)['xyz_bb']

    if not np.all(np.array([len(d['xyz_bb']) == len(xyz_ref) for d in datadicts], dtype=bool)):
        raise ValueError("All input PDB files must have the same number of CA atoms.")
    
    xyz_trg = np.stack([d['xyz_bb'] for d in datadicts])
    xyz_sup = superpose(xyz_ref, xyz_trg)

    for xyz, name, data in zip(xyz_sup, names, datadicts):
        data['xyz_bb'] = xyz
        data['xyz_ca'] = None
        new_file_name = args.output_dir + name + "_sup.pdb"
        with open(new_file_name, 'w') as f:
            f.write(write_pdb(**data))


if __name__ == '__main__':
    from pathlib import Path
    import sys
    sys.path.append(str(Path(__file__).parent.parent))

    #pdb_rmsd()
    pdb_superpose()