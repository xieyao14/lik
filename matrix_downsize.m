function B = matrix_downsize( A, downfactor)


assert( floor(downfactor)==downfactor) %assert downfactor to be integer

[m,n] = size(A);
assert( mod(m,downfactor)==0);
assert( mod(n,downfactor)==0);

m1 = m/downfactor;
B = reshape(A, [downfactor,m1,n]);
B = mean(B,1);
B = reshape(B, [m1,n]);


n1 = n/downfactor;
B = reshape(B', [downfactor,n1,m1]);
B = mean(B,1);
B = reshape(B, [n1,m1]);
B = B';

return;