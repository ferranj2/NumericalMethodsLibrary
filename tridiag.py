"""
%% Tridiagonal System solver v1.0
Written by J.A. Ferrand B.Sc (ID: 2431646)
Embry-Riddle Aeronautical University - Daytona Beach

%% Description
A simple function that solves linear systems of the form $Ax = b$ using
the Thomas algorithm. The function, therefore, is intended for linear
algebra problems where the matrix $A$ is tridiagonal. The $b$ vector is
referred to as the Right-Hand-Side (RHS). The $A$ matrix can be input in
one of two ways. Firstly as an $r \times r$ square matrix (which should
be avoided). Secondly (and preferrably) as an $r \times 3$ rectangular
matrix in which the first column is a zero that is vertically
concatenated with A's subdiagonal, the second column is the main
diagonal, and the third column is the superdiagonal vertically
concatenated with a zero. This 2nd method takes advantage of the inherent
sparsity of tridiagonal systems to save memory and improve speed. If the
entire matrix is input, then the code first checks to ensure that the
matrix is in fact tridiagonal before implementing the Thomas algorithm.

%% Changelog
v1.1,(05/15/2021): Improved index vetorization. Enabled multiple RHS.
v1.0,(05/13/2021): Initial Release (Single RHS only).

%% Syntax
* INPUT(*A*): Tridiagonal Coefficient Matrix of size $r$.
* INPUT(*b*): Right-Hand Side. (matrix that is $r\times n$)
* OUTPUT(*x*): Solution. (Same size as *b*)

%% Function Definition
"""
def tridiag(A,B):
    import numpy as np
    dA = A.shape
    rA,cA = dA[0],dA[1] #Get rows and columns of coefficient matrix.
    a = A.ravel() #Linearized reference for A matrix. (Vectorization purposes)
    dB = B.shape
    rB,cB = dB[0],dB[1] #Get rows and columns of Right-Hand-Side vector.
    b = B.ravel() #Linearized reference for RHS vector. (Vectorization purposes)
    if rA != rB:
        raise Exception('Coefficient matrix and RHS do not have equal rows!')
    if rA == cA: #Square Mode (rA by rA "A" matrix)
        print('Non-sparse matrices take longer to process than sparse ones!')
        nonzero = a != 0 #Count number of nonzero elements in the coefficient matrix.
        if sum(nonzero) > 3*rA-2: #If more zeroes than tridiagonal elements permit.
            raise Exception('Input matrix is NOT tridiagonal.')
        start = 0 #LI of first main diagonal element in square mode.
        dia = cA + 1 #LI for a main diagonal step in square mode.
        edpA = rA*cA-1 #LI of the last entry in the A matrix in square mode.
    elif rA != cA and cA == 3: #Sparse Mode (rA by 3 "A" matrix).
        start = 1 #LI of first main diagonal element in sparse mode.
        dia = cA #LI for a main diagonal step in sparse mode.
        edpA = cA*(rA-1)+1 #LI of the last entry in the A matrix in sparse mode.
    else:
        raise Exception('Invalid coefficient matrix. Must be square or a 3-column sparse matrix.')  
    
    #Thomas algorithm deployment
    X = np.zeros([rA,cB]) #Preallocate memory for the solution vector(s).
    x = X.ravel() #Linearized reference for solution vector(s).
    edpB = rB*cB-1 #LI of the last entry in the RHS.
    
    #Forward elimination of subdiagonal.
    row = np.arange(0,cB) #LI's of the first row (x & b)
    for i in range(0,rA-1):
        Mi = start + dia*i #LI of "ith" row superdiagonal element.
        Ui = Mi + 1 #LI of "ith+1" row main diagonal element.
        Mip1 = Mi + dia #LI of "ith" row main diagonal element.
        Lip1 = Mi + cA #LI of "ith+1" subdiagonal element of "ith" row.
        LoM = a[Lip1]/a[Mi] #Row elimination ratio.
        prow = row #Copy last row LI's
        row = prow + cB #Update LI's for next row.
        a[Mip1] = a[Mip1] - a[Ui]*LoM #Forward elimination update on main diagonal.
        b[row] = b[row] - b[prow]*LoM #Forward elimination on RHS vector(s).
    
    #Backward substitution to solve for x.
    row = np.arange(edpB-cB+1,edpB+1) #LI of the last row (x & b).
    x[row] = b[row]/a[edpA] #Solve for the last x-value.
    for i in range(1,rA):
        Mi = edpA - dia*i #LI of main diagonal element of "ith" row.
        Ui = Mi + 1 #LI of "ith" row superdiagonal element.
        prow = row[:] #Store indices of previous row.
        row = prow[:] - cB #Get indices of next row.
        x[row] = (b[row] - x[prow]*a[Ui])/a[Mi] #Solve for all "ith" x-values.      
    return X

#Testing the function
#import numpy as np
#A = np.array([[2,1,0,0],[1,2,1,0],[0,1,2,1],[0,0,1,2]], dtype = np.float64)
#B = np.array([[1,2],[1,2],[1,2],[1,2]], dtype = np.float64)
#X = tridiag(A,B)