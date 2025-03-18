# Newton-Raphson Algorithm
# f: scalar, gradient, Hessian
# x: initial guessing position
function [x, iter, xs] = newton(f, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g, H] = feval(f, x);
    if norm(g) < tol
      return
    endif
    # check degeneration
    if rank(H) < rows(x)
      return
    endif
    # searching direction
    s = -H\g;
    if norm(s) < tol
      break
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, max_iter);
    elseif
      a = 1.;
    endif
    # update position
    x += s*a;
  endfor

endfunction

