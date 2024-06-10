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
    from math import cos, atan2, atan
    if A.shape[0] != A.shape[1]:
        raise Exception("Matrix must be square!")
    
    
    N_points = 1_000
    Rs = off_diag_sum(A)
    theta = np.linspace(0, 2*np.pi, N_points)
    N = A.shape[0]
    
    # TODO: Optimize this shit here, disgusting disgusting for loops
    # This is going to be slow as fuck for large matrices
    # Potential ideas:
    #   Rewrite this as a library in CPP or some other 'fast' language
    #   Find another way to compute the boundary, (research alternative ways to do this)
    #   Potential for vectorization (Time to RTFM, fml)***
    #   What does industry do?
    for i in range(N):
        for j in range(i+1,N):
            """
            Alright, you just gotta believe me.  The below shit may appear like black magic sorcery and I 
            suspect is completely inaccessible to anyone else and will probably be inaccessible to me if I ever
            decide to come back to this.  So here goes trying to explain this sorcery...
            
            Future reference
            We are solving the following quartic polynomial with two dependent variables, $a and $b:
            r^4 - 2(a^2)(r^2)(cos 2 theta) - b^4 - a^4
            
            a and b are dependent variables, we can thus express one in terms of the other.
            
            We define new variable $alpha which makes our life easier somehow, but we gotta recover $a and $b at some point
            We do this through a series of steps, we solve a closely related matrix
            
            We obtain the related matrix by doing steps 1-3
            We get our alpha.
            
            Then, we gotta undo steps 3-1 to recover our a's and b's (incode: written as zeta and beta)
            """
            r1 = A[i][i] # Diagonal element; the unmodified complex form
            r1_ang = -atan(r1.imag / r1.real) # the angle form of r1
            if np.isnan(r1_ang):
                r1_ang = 0
            
            
            # Perform a shift
            A1 = A - r1*np.eye(N) # 1 - Shift, make A[i][i] zero
            
            r2 = -A1[j][j] # The unmodified complex form
            r2_ang = -atan(r2.imag / r2.real) # the angle of r2
            
            if np.isnan(r2_ang):
                r2_ang = 0
            
            A2 = r2*np.eye(N)*A1 # 2 - Rotate the diagonal A[j,j] entry towards the real line
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
                    z = r2*cos(t) + r2*sin(t)*j # Get the polar form
                    flag = z.real > 0
                    
                    # Intermediate step in undoing step 3 'unshift'
                    t = z + alpha # Intermediate value
                    beta = atan2(t.imag, t.real) 
                    
                    # Intermediate step in undoing 2 'rotation'
                    t = z - alpha
                    zeta = atan2(t.imag, t.real)
                    
                    # 
                    z = np.conj(r1_ang) * (np.conj(r2) * (z + alpha) + r1) 
                    
                    if alpha > 0: # Implies some sort of mirroring (flipped along some axis)
                        # In which case, we swap gamma and beta
                        beta, zeta = zeta, beta
                    
                        
                    return np.array([z, beta, zeta, flag])
                
                
    
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
    
    
    
    