%% LU Decomposition Algorithm
function[L, U, P]=luFactor(A)
% This funciton partially pivots and defactorizes a square matrix A into L 
% and U which are lower and upper triangular matrices. P is the permutation
% matrix, which keeps track of pivoting throuought the process.
% Inputs:
% A (square matrix)
% Outputs
% L (lower triangular matrix of A) 
% U (upper triangular matrix of A
% P (permutation matrix)
if nargin ~=1
    error('This function only requires one matrix input')
end
[m,n]=size(A);
if m~=n
    error('Input a square matrix');
end
P=eye(m); % Defining permutation matrix as an identity matrix
L=eye(m); % Defining lower triangular matrix as an identity matrix
U=A; % U is A prior to pivoting and Gaussian elimination

% Pivoting
    for i=1:m
       [value, row]=max(abs(U(i:m,i))); % 'row' is the row with the largest coefficient 
       row=row+i-1; % Determines if the row of largest coefficient is the same as the iteration
       if row~=i % Pivoting is required if row of largest coefficient does not equal iteration
           % Pivoting rows in U matrix
           u=U(i,:); % Temporary array to preserve rows in U matrix prior to pivoting
           U(i,:)=U(row,:);
           U(row,:)=u;
           % Pivoting rows in P matrix
           p=P(i,:); % Temporary array to preserve rows in P matrix prior to pivoting
           P(i,:)=P(row,:);
           P(row,:)=p;
           if i>=2 % L matrix requires no pivoting on first iteration
               % Pivoting matrix L
               l=L(i,1:i-1);
               L(i,1:i-1)=L(row,1:i-1);
               L(row,1:i-1)=l;
           end
       end
       % Gaussian elimination for L and U matrices
       for j=i+1:m % The first row is unchanged, so one less iteration is required
           L(j,i)=U(j,i)/U(i,i); % L matrix is defined by ratio U(j,i)/U(i,i)
           U(j,:)=U(j,:)-L(j,i)*U(i,:); % Eliminating variables in U matrix
       end
    end
    display(L)
    display(U)
    display(P)
end
