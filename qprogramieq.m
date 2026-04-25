# Quadratic programming by Active Set method
#   min |A*x - b|^2/2
#   s.t. C*x <= d
function [x, iter, xs] = qprogramieq(A, b, C, d, x, tol, max_iter)

  nv = rows(b); # number of variables
  nc = rows(d); # number of constraints
  ac = (C*x - d) >= 0; # active set
  xs = []; # searching history

  # With C_k*x_k = d_k, C_k*x_{k+1} = d_k and
  # x_{k+1} = x_k + s, the subproblem w.r.t
  # the search direction s_k is defined as:
  #   min |A*(x_k + s_k) - b|^2/2
  #   s.t. C_k*s_k = 0
  # where C_k, d_k are the current active set.

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # solve the subproblem of current active set
    id = find(ac); # ids of active constraints
    Ce = C(id, :);
    de = d(id, :);
    s = qprogrameq(A, b, Ce, de) - x;
    if norm(s) < tol
      return
    else # find the step not exceeding the set
      an = d - C*x; ad = C*s;
      a = ifelse(!ac & ad > 0, an./ad, 1);
      [a, k] = min(a);
      x += s*a;
      if a < 1 # there is a blocking constraint
        ac(k) = 1;
      endif
    endif
  endfor

endfunction

