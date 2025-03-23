# Quadratic Penalty method
function [x, iter, xs] = quadratic_penalty(f, c_eq, c_ieq, x, h, r, tol, max_iter)

  xs = []; # searching history

  function [u, g, H] = P(_x)
    [u, g, H] = penalty(f, c_eq, c_ieq, _x, h);
  endfunction

  x_p = x;

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g] = P(x);
    if norm(g) < tol
      return
    endif
    # solve unconstraint subproblem
    x = newton(@P, x, tol, max_iter);
    #x = bfgs(@P, x, tol, max_iter);
    if norm(x - x_p) < tol
      return
    endif
    # increase penalty factor
    h *= r;
    # backup results
    x_p = x;
  endfor

endfunction

#
function [u, g, H] = penalty(f, c_eq, c_ieq, x, h)

  n_f = nargout(f);

  if n_f >= 3
    [u, g, H] = feval(f, x);
  elseif n_f >= 2
    [u, g] = feval(f, x);
  else
    u = feval(f, x);
  endif

  if nargout(c_eq) >= 2
    [v, J] = feval(c_eq, x);
    u += v.^2*h*.5;
    if n_f >= 2
      g += J'*v*h;
    endif
    if n_f >= 3
      H += J'*J*h;
    endif
  elseif nargout(c_eq) >= 1
    v = feval(c_eq, x);
    u += v.^2*h*.5;
  endif

  if nargout(c_ieq) >= 2
    [v, J] = feval(c_ieq, x);
    v = max(v, 0);
    u += v.^2*h*.5;
    if n_f >= 2
      g += J'*v*h;
    endif
    if n_f >= 3
      s = sign(v);
      J.*= s*s';
      H += J'*J*h;
    endif
  elseif nargout(c_ieq) >= 1
    v = max(feval(c_ieq, x), 0);
    u += v.^2*h*.5;
  endif

endfunction

