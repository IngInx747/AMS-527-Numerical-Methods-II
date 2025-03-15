# Muller's method
function [x, iter] = solve_muller(f, x0, x1, x, eps, max_iter)

  # (f_k - f_{k-1}) / (x_k - k_{k-1})
  D = @(f, x1, x0)((f(x1) - f(x0))/(x1 - x0));

  for iter = 1 : max_iter
    h0 = x  - x1;
    h1 = x1 - x0;
    d0 = D(f, x , x1);
    d1 = D(f, x1, x0);
    a = (d0 - d1)/(h0 + h1);
    b = h0*a + d0;
    c = f(x);
    if norm(c) < eps
      break
    endif
    s0 = b + sqrt(b^2 - a*c*4);
    s1 = b - sqrt(b^2 - a*c*4);
    if norm(s0) > norm(s1)
      s = s0;
    else
      s = s1;
    endif
    s = -c*2/s;
    x0 = x1;
    x1 = x;
    x += s;
    if norm(s) < eps
      break
    endif
  endfor

endfunction
