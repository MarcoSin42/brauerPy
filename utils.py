import numpy as np

def off_diag_sum(A):
    """ Computes the sum of the off diagonal entries per row and returns it in the form of a vector

    Args:
        A (np.typing.ArrayLike): A (possibly complex) square matrix

    Returns:
        vector: Contains the off diagonal sums per row of a matrix
    """
    if A.shape[0] != A.shape[1]:
        raise Exception("Matrix must be square!")
    
    A_mod = abs(A)
    
    return A_mod.sum(axis=1) - np.diag(A_mod)

def eigs(A):
    """ Computes the eigenvalues of a (possibly complex-valued) matrix and returns real and imaginary part into two vectors

    Args:
        A (np.typing.ArrayLike): (possibly complex-valued) matrix
    
    Returns:
        2xn array: first row contains the real part and second row contains the imaginary part
    """
    
    if A.shape[0] != A.shape[1]:
        raise Exception("Matrix must be square!")
    
    e = np.linalg.eigvals(A)
    return e.real, e.imag

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    print("Running test")
    
    a = np.vander([1,2,3], 3)
    
    cx = a.diagonal().real
    cy = a.diagonal().imag
    
    radii = off_diag_sum(a)
    n = cx.shape[0]
    
    fig, ax = plt.subplots()
    
    bbox = (-20, 20)
    ax.set_xlim(bbox)
    ax.set_ylim(bbox)
    
    for i in range(n):
        ax.add_patch(
            plt.Circle(
                (cx[i], cy[i]),
                radii[i],
                color=(0,0,1,0.1)
            )
    )
    
    r, i = eigs(a)
    ax.plot(r,i, 'ro')
    
    
    
    
    fig.savefig("3x3_testimg")
    
    
    
    