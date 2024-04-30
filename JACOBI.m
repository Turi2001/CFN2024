function VEC = JACOBI(xhat,A,b)
  VEC = [];
  VEC = [VEC xhat];
  index = 0;
  D = diag(diag(A));
  L = tril(A, -1);
  U = triu(A, 1);
  for index = 1:20;
     xhat = inv(D)*(b-(L-U)*xhat);
     VEC = [VEC xhat];
  endfor
endfunction