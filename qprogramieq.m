# Quadratic programming by Active Set method
#   min x'*A*x/2 + x'*b
#   s.t. C*x <= d
function [x, iter, xs] = qprogramieq(A, b, C, d, x, tol, max_iter)

  nv = rows(b); # number of variables
  nc = rows(d); # number of constraints
  ac = (C*x - d) >= 0; # active set bits
  xs = []; # searching history

  # The subproblem regarding the current active set:
  #   min x_k'*A*x_k/2 + x_k'*b
  #   s.t. C_k * x_k = d_k
  # where {k \in K} are ids of the active equalities.

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # solve the subproblem of current active set
    id = find(ac); # ids of active constraints
    Ce = C(id, :);
    de = d(id, :);
    [s, y] = qprogrameq(A, b, Ce, de); s -= x;
    if norm(s) < tol
      if isempty(y) || min(y) > -tol # reach a critical point
        return
      else # improve by dropping one of the active constraints
        [_, k] = min(y);
        ac(id(k)) = 0;
      endif
    else # find the step not breaking the inactive constraints
      an = d - C*x; ad = C*s;
      a = ifelse(!ac & ad > 0, an./ad, 1);
      [a, k] = min(a);
      x += s*a;
      if a < 1 # add a blocking constraint to the active set
        ac(k) = 1;
      endif
    endif
  endfor

endfunction

