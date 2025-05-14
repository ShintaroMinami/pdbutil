import numpy as np

def kabsch(A, B):
    """
    Kabsch algorithm to find the optimal rotation matrix that minimizes the RMSD
    between two sets of points A and B.
    Parameters:
    A : numpy array of shape (N, 3) or (1, N, 3)
        First set of points.
    B : numpy array of shape (N, 3) or (1, N, 3)
        Second set of points.
    Returns:
    numpy array of shape (N, 3)
        Rotated version of A that best matches B.
    """
    Y = np.expand_dims(A, axis=0) if len(A.shape) == 2 else A
    X = np.expand_dims(B, axis=0) if len(B.shape) == 2 else B
    X = X - X.mean(axis=-2, keepdims=True)
    Y_mean = Y.mean(axis=-2, keepdims=True)
    Y = Y - Y_mean
    C = np.einsum('b i j, b j k -> b i k', X.transpose(0,2,1), Y)
    V, _, W = np.linalg.svd(C)
    # det sign for direction correction
    d = np.sign(np.linalg.det(V) * np.linalg.det(W))
    V[:,:,-1] = V[:,:,-1] * np.repeat(np.expand_dims(d, axis=1), 3, axis=1)
    # calc rotation
    R = np.einsum('b i j, b j k -> b i k', V, W)
    # rotate
    X = np.einsum('b a c, b c r -> b a r', X, R)
    return (Y+Y_mean).squeeze(), (X+Y_mean).squeeze()


def kabsch2(A, B):
    """
    Kabsch algorithm to find the optimal rotation matrix that minimizes the RMSD
    between two sets of points A and B.
    Parameters:
    A : numpy array of shape (N, 3) or (1, N, 3)
        First set of points.
    B : numpy array of shape (N, 3) or (1, N, 3)
        Second set of points.
    Returns:
    numpy array of shape (N, 3)
        Rotated version of B that best matches A.
    """
    X = np.expand_dims(A, axis=0) if len(A.shape) == 2 else A
    Y = np.expand_dims(B, axis=0) if len(B.shape) == 2 else B
    Y = Y - Y.mean(axis=-2, keepdims=True)
    X_mean = X.mean(axis=-2, keepdims=True)
    X = X - X_mean
    C = np.einsum('b i j, b j k -> b i k', Y.transpose(0,2,1), X)
    V, _, W = np.linalg.svd(C)
    # det sign for direction correction
    d = np.sign(np.linalg.det(V) * np.linalg.det(W))
    V[:,:,-1] = V[:,:,-1] * np.repeat(np.expand_dims(d, axis=1), 3, axis=1)
    # calc rotation
    R = np.einsum('b i j, b j k -> b i k', V, W)
    # rotate
    Y = np.einsum('b a c, b c r -> b a r', Y, R)
    return (X+X_mean).squeeze(), (Y+X_mean).squeeze()


def rmsd(A, B, return_xyz=False):
    X, Y = kabsch(A, B)
    rmsd = np.sqrt(((X-Y)**2).sum(axis=-1).mean(axis=-1))
    return (rmsd, X, Y) if return_xyz else rmsd

