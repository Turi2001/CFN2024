function x = JACOBI(A,b)
  % ATTENTION
  % La fonction ne permet de resoudre le système que si la matrice associé 
  % est diagonal dominante 
   n = length(b);
   x = zeros(n,1);
   
   L = tril(A, -1);
   U = triu(A, 1);
   D = diag(A);
   
    max_iter = 100;
    tol = 1e-6;
    iter = 0;
    
   while iter < max_iter;
     x_new = (b - (L + U) *x) ./ D;
     if norm(x_new-x) < tol
       x = x_new;
       break;
     endif
     x = x_new;
     iter = iter +1;
   endwhile
   if iter == max_iter
        warning('La méthode de Jacobi n''a pas convergé dans le nombre maximal d''itérations.');
   end
   
endfunction
