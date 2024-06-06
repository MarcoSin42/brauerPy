import numpy as np
from math import sin,cos

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

def brauer(A):
    from math import cos
    if A.shape[0] != A.shape[1]:
        raise Exception("Matrix must be square!")
    
    
    N_points = 1_000
    Rs = off_diag_sum(A)
    theta = np.linspace(0, 2*np.pi, N_points)
    N = A.shape[0]
    
    # TODO: Optimize this shit here
    # This is going to be slow as fuck for large matrices
    # Potential ideas:
    #   Rewrite this as a library in CPP or some other 'fast' language
    #   Find another way to compute the boundary, (research alternative ways to do this)
    #   Potential for vectorization (Time to RTFM, fml)***
    #   What does industry do?
    for i in range(N):
        for j in range(i+1,N):
            print("ij")
            w = A[i][i]
            
            # Perform a shift
            A1 = A - w*np.eye(N) # 1 - Shift
            gamma = -A1[j][j]
            
            A2 = gamma*np.eye(N)*A1 # 2 - Rotate the diagonal A[j,j] entry towards the real line
            alpha = A2[j][j]/2
            A3 = A2 - alpha*np.eye(N) # 3 - Shift again
            
            for t in theta:
                P = np.polynomial.Polynomial((
                    alpha**4 - (Rs[i]*Rs[j])**2, # Constant coef
                    0, # x
                    -2*alpha**2*cos(2*t), # x^2
                    0, # x^3
                    1  # x^4
                ))
                
                tol = 0.0001 
                soln = P.roots() # Solutions to our rotated and shifted matrix
                
                real_roots = soln[soln.imag < tol] # Retrieve real roots only! 
                
                for r in real_roots:
                    z = w*cos(t) + w*sin(t)*j
                    flag = z.real > 0
                    
                    
                
                
    
    return None
    
    

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
    
    brauer(a)
    
    
    
    
    fig.savefig("3x3_testimg")
    
    
    
    