%% Tridiagonal System solver v1.0
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%% Description
% A simple function that solves linear systems of the form $Ax = b$ using 
% the Thomas algorithm. The function, therefore, is intended for linear 
% algebra problems where the matrix $A$ is tridiagonal. The $b$ vector is 
% referred to as the Right-Hand-Side (RHS). The $A$ matrix can be input in 
% one of two ways. Firstly as an $r \times r$ square matrix (which should
% be avoided). Secondly (and preferrably) as an $r \times 3$ rectangular
% matrix in which the first column is a zero that is vertically 
% concatenated with A's subdiagonal, the second column is the main 
% diagonal, and the third column is the superdiagonal vertically 
% concatenated with a zero. This 2nd method takes advantage of the inherent
% sparsity of tridiagonal systems to save memory and improve speed. If the 
% entire matrix is input, then the code first checks to ensure that the 
% matrix is in fact tridiagonal before implementing the Thomas algorithm.
%%
% Sample Square ("entire") tridiagonal matrix
%%
% $$A = \left[\begin{array}{cccc}2&1&0&0\\1&2&1&0\\0&1&2&1\\0&0&1&2\end{array}\right]$
%%
% Sample Sparse matrix format (zeros are placeholder values)
%%
% $$ A = \left[\begin{array}{ccc}0&2&1\\1&2&1\\1&2&1\\1&2&0\end{array}\right]$
%% Changelog
%  v1.1,(05/15/2021): Improved index vetorization. Enabled multiple RHS.
%  v1.0,(05/13/2021): Initial Release (Single RHS only).
%% Syntax
% * INPUT(*A*): Tridiagonal Coefficient Matrix of size $r$. 
% * INPUT(*b*): Right-Hand Side. (matrix that is $r\times n$) 
% * OUTPUT(*x*): Solution. (Same size as *b*)
%% Function Definition
function x = tridiag(A,B)
[rA,cA] = size(A); %Get rows and columns of coefficient matrix.
[rB,cB] = size(B); %Get rows and columns of Right-Hand-Side vector.
if rA ~= rB
    error('Mismatch between matrix and RHS vector rows.')
elseif cB > 1
    warning('Multiple RHS detected!')
end
x = zeros(rA,cB); %Preallocate memory for the solution vector.
dia = rA + 1; %Linear index for a main diagonal step.
edpA = rA*cA; %Linear index of the last entry in the A matrix.
edpB = rB*cB; %Linear index of the last entry in the RHS.
trb = edpB - rB + 1; %Linear index of top-right element in b.
if rA == cA %Square matrix supplied
    warning('Non-sparse matrices take longer to process than sparse ones!')
    nonzero = A ~= 0; %Count number of nonzero elements in the coefficient matrix.
    if sum(nonzero(:)) > 3*rA-2 %If more zeroes than tridiagonal elements permit.
        error('Input matrix is NOT tridiagonal.')
    else %Deploy the thomas algorithm.
        row = 1:rA:trb; %Linear indices of the first row (x & b)
        for i = 1:rA-1 %Forward elimination of subdiagonal.
            idia = i*dia; %Diagonal step times "i"
            idiap1 = idia + 1; %"idia" plus 1.
            M = idiap1-dia; %Linear index of the main diagonal element of the "ith" column.
            L = M+1; %Linear index of the subdiagonal element of the "ith" column.
            LoM = A(L)/A(M); %Row elimination ratio.
            prow = row;
            row = prow + 1;
            B(row) = B(row) - B(prow)*LoM; %Forward elimination on RHS vector(s).
            A(idiap1) = A(idiap1) - A(idia)*LoM; %Forward elimination update on main diagonal.
        end
        %Backward substitution to solve for x.
        row = rA:rA:edpB; %Linear indices of the last row (x & b).
        x(row) = B(row)/A(edpA); %Solve for the last x-value.
        for i = 1:rA-1 
            prow = row; %Store indices of previous row.
            row = row - 1; %Get indices of next row.
            M = edpA - dia*i; %Linear index of main diagonal element of "ith" row.
            U = M + rA; %Linear index of superdiagonal element.
            x(row) = (B(row) - x(prow)*A(U))/A(M); %Solve for all "ith" x-values.
        end
    end
elseif rA~=cA && cA ~= 3 %Invalid matrix input.
    error('Invalid Input matrix. Must be square or a 3-column matrix.')
else %Tridiagonal system input as the r by 3 matrix. (PREFERRED)
    row = 1:rA:trb; %First row linear indices (x & b).
    for i = 2:rA
        prow = row; %Store previous row indices.
        row = prow + 1; %Update row indices.
        M = i+rA; %Linear index of the main diagonal element.
        L = M-1; %Linear index of the subdiagonal element.
        A(M) = A(M) - A(i)*A(L+rA)/A(L); 
        B(row) = B(row) - A(i)*B(prow)/A(L);
    end
    trA=2*rA; %"two times rA"
    row = rA:rA:edpB; %Last row linear indices (x & b).
    x(row) = B(row)/A(trA); %Compute the very last x-value.
    for i = rA-1:-1:1
        prow = row; %Store previous row indices.
        row = prow -1; %Update row indices.
        x(row) = (B(row) - A(i+trA)*x(prow))/A(i+rA); %Use backward substitution to solve remaining x.
    end
end
end
%% Sources
% * "Numerical Methods for Engineers and Scientists" 2nd ed. (2001) by Joe
% D. Hoffman
% * "Essential Computational Fluid Dynamics" 2nd ed. (2010) by Oleg Zikanov