# Broyden's Quasi-Newton method (direct Hessian)
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
function [x, iter, xs] = broydend(f, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  if nargout(f) >= 3
    [~, g, B] = feval(f, x); # initial gradient and Hessian
  elseif nargout(f) >= 2
    [~, g] = feval(f, x); # initial gradient
    B = eye(rows(x));
  else
    error("Target function requires C1 smoothness!");
  endif

  if rank(B) < rows(x)
    error("Initial Hessian is degenerated!");
  endif

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # searching direction
    s = -B\g;
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
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      return
    endif
    # update Hessian surrogate
    B += g*s' / (s'*s);
    if rank(B) < rows(x)
      error("Hessian is degenerated!");
    endif
  endfor

endfunction

