# Muller's method
# f: target function/equation
# [x, x0, x1]: initial guess/interval
function [x, iter, xs] = muller(f, x, x0, x1, tol, max_iter)

  xs = [x0, x1]; # searching history

  f0 = feval(f, x0);
  f1 = feval(f, x1);

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    fx = feval(f, x);
    # interpolation
    h0 = x1 - x0;
    h1 = x  - x1;
    d0 = (f1 - f0)/h0;
    d1 = (fx - f1)/h1;
    a = (d1 - d0)/(h0 + h1);
    b = h1*a + d1;
    c = fx;
    if norm(c) < tol
      return
    endif
    # searching direction
    s0 = b + sqrt(b^2 - a*c*4);
    s1 = b - sqrt(b^2 - a*c*4);
    if norm(s0) > norm(s1)
      s = s0;
    else
      s = s1;
    endif
    s = -c*2/s;
    if norm(s) < tol
      return
    endif
    # backup results
    x0 = x1;
    x1 = x;
    f0 = f1;
    f1 = fx;
    # update position
    x += s;
  endfor

endfunction

