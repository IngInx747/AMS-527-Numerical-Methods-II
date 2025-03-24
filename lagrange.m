# Augmented Lagrangian method
function [x, iter, xs] = lagrange(f, c_eq, c_ieq, x, h, r, tol, max_iter)

endfunction

# penalty of inequality constraints
function [u, g, H] = penalty(c, x)

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

