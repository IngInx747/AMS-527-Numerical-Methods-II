# Augmented Lagrangian method
# f: target function
# ce: equality constraints
# ci: inequality constraints
function [x, iter, xs] = lagrange(f, ce, ci, x, r, tol, max_iter)

  xs = []; # searching history

  # Lagrange multipliers
  if nargout(ce) > 0
    ze = zeros(size(ce(x)));
  else
    ze = [];
  endif

  if nargout(ci) > 0
    zi = zeros(size(ci(x)));
  else
    zi = [];
  endif

  # define subproblem
  function [u, g, H] = sub(_x)
    [u, g, H] = subproblem(f, ce, ci, _x, ze, zi, r);
  endfunction

  x_p = x;

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g] = sub(x);
    if norm(g) < tol
      return
    endif
    # solve unconstraint subproblem
    #x = newton(@sub, x, tol, max_iter);
    x = bfgs(@sub, x, tol, max_iter);
    if norm(x - x_p) < tol
      return
    endif
    # update multipliers
    if nargout(ce) > 0
      ze += ce(x)*r;
    endif
    if nargout(ci) > 0
      zi = max(zi + ci(x)*r, 0);
    endif
    # backup results
    x_p = x;
  endfor

endfunction

# target + dual + penalty
function [u, g, H] = subproblem(f, ce, ci, x, ze, zi, r)

  n_f = nargout(f);

  if n_f >= 3
    [u, g, H] = feval(f, x);
  elseif n_f >= 2
    [u, g] = feval(f, x);
  else
    u = feval(f, x);
  endif

  if nargout(ce) >= 2
    [v, J] = feval(ce, x);
    # multiplier terms
    u += v'*ze;
    if n_f >= 2
      g += J'*ze;
    endif
    # penalty terms
    u += v'*v*r*.5;
    if n_f >= 2
      g += J'*v*r;
    endif
    if n_f >= 3
      H += J'*J*r;
    endif
  elseif nargout(ce) >= 1
    v = feval(ce, x);
    # multiplier terms
    u += v'*ze;
    # penalty terms
    u += v'*v*r*.5;
  endif

  if nargout(ci) >= 2
    [v, J] = feval(ci, x);
    # multiplier terms
    v = v + zi/r;
    u -= zi'*zi/r*.5;
    # masking
    v = max(v, 0);
    J .*= sign(v);
    # penalty terms
    u += v'*v*r*.5;
    if n_f >= 2
      g += J'*v*r;
    endif
    if n_f >= 3
      H += J'*J*r;
    endif
  elseif nargout(ci) >= 1
    v = feval(ci, x);
    # multiplier terms
    v = v + zi/r;
    u -= zi'*zi/r*.5;
    # masking
    v = max(v, 0);
    # penalty terms
    u += v'*v*r*.5;
  endif

endfunction

