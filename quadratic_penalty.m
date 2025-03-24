# Quadratic Penalty method
function [x, iter, xs] = quadratic_penalty(f, c_eq, c_ieq, x, h, r, tol, max_iter)

  xs = []; # searching history

  function [u, g, H] = P(_x)
    [u, g, H] = subproblem(f, c_eq, c_ieq, _x, h);
  endfunction

  function [u, g, H] = P_eq(_x)
    u = penalty_eq(c_eq, _x);
  endfunction

  function [u, g] = P_ieq(_x)
    u = penalty_ieq(c_ieq, _x);
  endfunction

  x_p = x;

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g] = P(x);
    p = (P_eq(x) + P_ieq(x))*h;
    if norm(g) < tol && p < tol
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

# target + penalty
function [u, g, H] = subproblem(f, c_eq, c_ieq, x, h)

  n_f = nargout(f);

  if n_f >= 3
    [u, g, H] = feval(f, x);
  elseif n_f >= 2
    [u, g] = feval(f, x);
  else
    u = feval(f, x);
  endif

  if nargout(c_eq) >= 2
    [c, d, B] = penalty_eq(c_eq, x);
    u += c*h;
    if n_f >= 2
      g += d*h;
    endif
    if n_f >= 3
      H += B*h;
    endif
  elseif nargout(c_eq) >= 1
    c = penalty_eq(c_eq, x);
    u += c*h;
  endif

  if nargout(c_ieq) >= 2
    [c, d, B] = penalty_ieq(c_ieq, x);
    u += c*h;
    if n_f >= 2
      g += d*h;
    endif
    if n_f >= 3
      H += B*h;
    endif
  elseif nargout(c_ieq) >= 1
    c = penalty_ieq(c_ieq, x);
    u += c*h;
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
    s = sign(v);
    J.*= s*s';
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

