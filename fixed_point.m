# Fixed-point iteration
function [x, iter, xs] = fixed_point(f, x, tol, max_iter)

  xs = []; # searching history

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # backup results
    x_0 = x;
    # update position
    x = feval(f, x);
    if norm(x - x_0) < tol
      break
    endif
  endfor

endfunction

