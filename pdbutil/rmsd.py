import numpy as np
from .pdbio import BB_ATOMS_TO_INDEX


def _rotate(X:np.array, rot:np.array, vec:np.array)->np.array:
    return np.einsum('b a c, b c r -> b a r', X, rot) + vec


def kabsch(A:np.array, B:np.array, return_rotmat=False):
    """
    Kabsch algorithm to find the optimal rotation matrix that minimizes the RMSD
    between two sets of points A and B.
    Parameters:
    A : numpy array of shape (N, 3) or (B, N, 3)
        First set of points.
    B : numpy array of shape (N, 3) or (B, N, 3)
        Second set of points.
    Returns:
    numpy array of shape (N, 3) or (B, N, 3)
        Rotated version of B that best matches A.
    """
    X = np.expand_dims(A, axis=0) if len(A.shape) == 2 else A
    Y = np.expand_dims(B, axis=0) if len(B.shape) == 2 else B
    X_mean = X.mean(axis=-2, keepdims=True)
    Y_mean = Y.mean(axis=-2, keepdims=True)
    C = np.einsum('b i j, b j k -> b i k', (Y-Y_mean).transpose(0,2,1), (X-X_mean))
    V, _, W = np.linalg.svd(C)
    # det sign for direction correction
    d = np.sign(np.linalg.det(V) * np.linalg.det(W))
    V[:,:,-1] = V[:,:,-1] * np.repeat(np.expand_dims(d, axis=1), 3, axis=1)
    # calc rotation
    rot = np.einsum('b i j, b j k -> b i k', V, W)
    vec = np.einsum('b a c, b c r -> b a r', -Y_mean, rot) + X_mean
    if return_rotmat:
        return rot, vec
    else:
        return X, _rotate(Y, rot, vec)


def _format_xyz(xyz_ref, xyz_trg, axis_L=None):
    xyz_ref = np.expand_dims(xyz_ref, axis=0) if xyz_ref.shape[0]!=1 else xyz_ref
    length = xyz_ref.shape[1]
    assert length > 4, "xyz_ref must have more than 4 residues"
    if xyz_ref.ndim == 3:
        xyz_ref_ca = xyz_ref
    elif xyz_ref.ndim == 4:
        xyz_ref_ca = xyz_ref[:,:,BB_ATOMS_TO_INDEX['CA']] # Only Ca atom
    else:
        raise ValueError("xyz_ref shape must be [L, 3], [1, L, 3], [L, 4, 3], or [1, L, 4, 3]")
    if axis_L is None:
        axis_L = np.where(np.array(xyz_trg.shape) == length)[0]
        assert len(axis_L) > 0, "Different length of xyz_reference and xyz_targets"
    if axis_L == 0:
        xyz_trg = np.expand_dims(xyz_trg, axis=0)
    if xyz_trg.ndim == 3:
        xyz_trg_ca = xyz_trg
    elif xyz_trg.ndim == 4:
        xyz_trg_ca = xyz_trg[...,1,:]
    else:
        raise ValueError("xyz_trg shape must be [L, 3], [B, L, 3], [L, 4, 3], or [B, L, 4, 3]")
    return xyz_trg, xyz_ref_ca, xyz_trg_ca


def superpose(
        xyz_reference:np.array,
        xyz_targets:np.array,
        return_rmsd:bool=False,
        axis_L:int=None
        ):
    """
    Superpose the target coordinates to the reference coordinates using Kabsch algorithm.
    Parameters:
    xyz_reference : numpy array of shape (L, 3) or (1, L, 3)
        Reference coordinates.
    xyz_targets : numpy array of shape (L, 3) or (B, L, 3)
        Target coordinates.
    return_rmsd : bool, optional
        If True, return the RMSD value. Default is False.
    axis_L : int, optional
        Axis of the length of the coordinates. Default is None.
    Returns:
    numpy array of shape (L, 3) or (B, L, 3)
        Superposed coordinates.
    float
        RMSD value if return_rmsd is True.
    """
    xyz_trg_format, xyz_ref_ca, xyz_trg_ca = _format_xyz(xyz_reference, xyz_targets, axis_L=axis_L)
    shape_org = xyz_targets.shape
    # Kabsch superposition
    rot, vec = kabsch(xyz_ref_ca, xyz_trg_ca, return_rotmat=True)
    # Flatten xyz
    xyz_trg_flat = xyz_trg_format.reshape(xyz_trg_format.shape[0], -1, 3)
    # Rotate and translate
    xyz_trg_flat_rot = _rotate(xyz_trg_flat, rot, vec)
    xyz_superposed = xyz_trg_flat_rot.reshape(shape_org)
    if return_rmsd:
        rmsd = np.sqrt(((xyz_ref_ca-_rotate(xyz_trg_ca, rot, vec))**2).sum(axis=-1).mean(axis=-1)).squeeze()
        return xyz_superposed, rmsd
    else:
        return xyz_superposed


def calc_rmsd(
        A:np.array,
        B:np.array
        )->np.array:
    """
    Calculate the RMSD between two sets of coordinates A and B.
    Parameters:
    A : numpy array of shape (L, 3) or (B1, L, 3)
        First set of coordinates.
    B : numpy array of shape (L, 3) or (B2, L, 3)
        Second set of coordinates.
    Returns:
    numpy array of shape (B1, B2)
        RMSD values between each pair of coordinates.
    """
    assert (A.ndim <= 3) and (B.ndim <= 3), "A and B must be 2D (L,3) or 3D (B,L,3) arrays of C-alpha coordinates"
    A = np.expand_dims(A, axis=0) if A.ndim == 2 else A
    B = np.expand_dims(B, axis=0) if B.ndim == 2 else B
    out_shape = (A.shape[0], B.shape[0])
    A_repeat = A[:,None].repeat(B.shape[0], axis=1).reshape(-1, A.shape[-2], 3)
    B_repeat = B[None,:].repeat(A.shape[0], axis=0).reshape(-1, B.shape[-2], 3)
    X, Y = kabsch(A_repeat, B_repeat)
    rmsd = np.sqrt(((X-Y)**2).sum(axis=-1).mean(axis=-1)).reshape(out_shape)
    return rmsd.squeeze()



