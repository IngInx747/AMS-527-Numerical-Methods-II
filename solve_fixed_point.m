#
function [x, iter] = solve_fixed_point(f, x, tol, max_iter)

  for iter = 1 : max_iter
    x_0 = x;
    x = f(x);
    if norm(x - x_0) < tol
      break
    endif
  endfor

endfunction

