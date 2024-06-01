import numpy as np

def off_diag_sum(A: np.typing.ArrayLike):
    """ Computes the sum of the off diagonal entries per row and returns it in the form of a vector

    Args:
        A (np.typing.ArrayLike): A (possibly complex) square matrix

    Returns:
        vector: Contains the off diagonal sums per row of a matrix
    """
    if A.shape[0] != A.shape[1] or len(A.shape):
        return
    
    A_mod = abs(A)
    
    return A_mod.sum(axis=1) - np.diag(A_mod)

def eigs(A: np.typing.ArrayLike):
    """ Computes the eigenvalues of a (possibly complex-valued) matrix and returns real and imaginary part into two vectors

    Args:
        A (np.typing.ArrayLike): (possibly complex-valued) matrix
    
    Returns:
        2xn array: first row contains the real part and second row contains the imaginary part
    """
    
    if A.shape[0] != A.shape[1] or len(A.shape):
        return
    
    return a.real, a.imaginary

if __name__ == '__main__':
    import matplotlib as plt
    
    print("Running test")
    
    