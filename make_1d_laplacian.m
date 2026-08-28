function L = make_1d_laplacian(n)

%return 1d lapalcian matrix L = D-A

A = spdiags( ones(n,1), [1], n,n);
A = A+A';

dA =sum(A,2);
L = diag(dA)-A;

end