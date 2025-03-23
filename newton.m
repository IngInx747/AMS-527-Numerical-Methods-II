# Newton-Raphson Algorithm
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
function [x, iter, xs] = newton(f, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  if nargout(f) < 3
    error("Target function requires C2 smoothness!");
  endif

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
      error("Hessian is degenerated!");
    endif
    # searching direction
    s = -H\g;
    if norm(s) < tol
      return
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, :);
    else
      a = 1.;
    endif
    # update position
    x += s*a;
  endfor

endfunction

