# Broyden–Fletcher–Goldfarb–Shanno method (inversed Hessian)
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
function [x, iter, xs] = bfgs(f, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  if nargout(f) >= 3
    [~, g, B] = feval(f, x); # initial gradient and Hessian
  elseif nargout(f) >= 2
    [~, g] = feval(f, x); # initial gradient
    B = eye(rows(x));
  else
    error("Target function requires C1 smoothness!");
  endif

  if !issymmetric(B)
    error("Initial Hessian is not symmetric!");
  endif

  C = inv(B); # inversed Hessian surrogate

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # searching direction
    s = -C*g;
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
    # backup results
    y = g;
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      return
    endif
    # delta gradient
    y = g - y;
    # update Hessian surrogate
    R = y*s'/(y'*s);
    C += R'*C*R - R'*C - C*R + s*s'/(y'*s);
    if rank(C) < rows(x)
      error("Hessian is degenerated!");
    endif
  endfor

endfunction

