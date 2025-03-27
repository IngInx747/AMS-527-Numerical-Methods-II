# Quadratic Penalty method
# f: target function
# ce: equality constraints
# ci: inequality constraints
# x: initial guess
function [x, iter, xs] = qpenalty(f, ce, ci, x, r, t, tol, max_iter)

  xs = []; # searching history

  function [u, g, H] = sub(_x)
    [u, g, H] = subproblem(f, ce, ci, _x, r);
  endfunction

  function u = pen(_x)
    u = penalty_eq(ce, _x) + penalty_ieq(ci, _x);
  endfunction

  x_p = x;

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g] = sub(x);
    if norm(g) < tol && pen(x)*r < tol
      return
    endif
    # solve unconstraint subproblem
    #x = newton(@sub, x, tol, max_iter);
    x = bfgs(@sub, x, tol, max_iter);
    if norm(x - x_p) < tol
      return
    endif
    # increase penalty factor
    r *= t;
    # backup results
    x_p = x;
  endfor

endfunction

# target + penalty
function [u, g, H] = subproblem(f, ce, ci, x, r)

  n_f = nargout(f);

  if n_f >= 3
    [u, g, H] = feval(f, x);
  elseif n_f >= 2
    [u, g] = feval(f, x);
  else
    u = feval(f, x);
  endif

  if nargout(ce) >= 2
    [c, d, B] = penalty_eq(ce, x);
    u += c*r;
    if n_f >= 2
      g += d*r;
    endif
    if n_f >= 3
      H += B*r;
    endif
  elseif nargout(ce) >= 1
    c = penalty_eq(ce, x);
    u += c*r;
  endif

  if nargout(ci) >= 2
    [c, d, B] = penalty_ieq(ci, x);
    u += c*r;
    if n_f >= 2
      g += d*r;
    endif
    if n_f >= 3
      H += B*r;
    endif
  elseif nargout(ci) >= 1
    c = penalty_ieq(ci, x);
    u += c*r;
  endif

endfunction

# penalty of equality constraints
function [u, g, H] = penalty_eq(c, x)

  if nargout(c) >= 2
    [v, J] = feval(c, x);
    u = v'*v*.5;
    g = J'*v;
    H = J'*J;
  elseif nargout(c) >= 1
    v = feval(c, x);
    u = v'*v*.5;
  else
    u = 0;
  endif

endfunction

# penalty of inequality constraints
function [u, g, H] = penalty_ieq(c, x)

  if nargout(c) >= 2
    [v, J] = feval(c, x);
    v = max(v, 0);
    J.*= sign(v);
    u = v'*v*.5;
    g = J'*v;
    H = J'*J;
  elseif nargout(c) >= 1
    v = feval(c, x);
    v = max(v, 0);
    u = v'*v*.5;
  else
    u = 0;
  endif

endfunction

