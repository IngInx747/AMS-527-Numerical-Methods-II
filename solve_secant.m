# Secant method
function [x, iter] = solve_secant(f, x0, x, tol, max_iter)

  # (f_k - f_{k-1}) / (x_k - k_{k-1})
  D = @(f, x1, x0)((f(x1) - f(x0))/(x1 - x0));

  for iter = 1 : max_iter
    y = f(x);
    if abs(y) < tol
      break
    endif
    d = D(f, x, x0);
    s = -y / d; # J(x_k)*s_k = -f(x_k)
    x0 = x;
    x += s; # s_k = x_{k+1} - x_k
    if abs(s) < tol
      break
    endif
  endfor

endfunction
