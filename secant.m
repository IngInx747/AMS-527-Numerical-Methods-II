# Secant method
# f: [scalar, gradient] (optimal problem)
#    [merit, equations] (solving equations)
# [x, x0]: initial guess/interval
function [x, iter, xs] = secant(f, x, x0, tol, max_iter)

  xs = [x0]; # searching history

  [~, g0] = feval(f, x0);

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      break
    endif
    # secant equation
    d = (g - g0)/(x - x0);
    # searching direction
    s = -g/d;
    if abs(s) < tol
      break
    endif
    # backup results
    x0 = x;
    g0 = g;
    # update position
    x += s;
  endfor

endfunction

