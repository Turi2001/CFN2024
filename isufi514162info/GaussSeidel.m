function x = GaussSeidel(A,b)
  n = length(b);
  x = zeros(n,1);
  
  L = tril(A,-1);
  U = triu(A,1);
  D = diag(diag(A));
  
  max_iter =100;
  tol = 1e-6;
  iter = 0;
  
  while iter<max_iter
    x_new = inv(D + L)*(b-U*x);
    if norm(x_new-x) < tol
       x = x_new;
       break;
     endif
     x = x_new;
     iter = iter+1;
  endwhile
  
  if iter == max_iter
        warning('La méthode de Jacobi n''a pas convergé dans le nombre maximal d''itérations.');
  end
endfunction
