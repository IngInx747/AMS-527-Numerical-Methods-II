# Quadratic programming
#   min |A*x - b|^2/2
#   s.t. C*x = d
function [x, v] = qprogrameq(A, b, C, d)

  x = zeros(size(b));

  # Reduce equality constraints
  #   C*x = d
  # by the transform
  #   x = T*z + q
  [T, q] = null_space(C, d, x);

  # The quadratic optimal problem becomes
  #   min  |A_*z - b_|^2/2
  #   where A_ = A*T
  #         b_ = b - A*q
  b -= A*q;
  A = A*T;

  # Solve the reduced problem
  z = (A'*A) \ (A'*b);
  x = T*z + q;

  # Lagrange multipliers (C*T is full rank)
  if nargout > 2
    v = (C*T)'\(T'*(A'*b - A'*A*x));
  endif

  # This method gives the same result as solved
  # by an argumented Lagrangian system
  # | A'A  C' | * | x | = | A'b |
  # |  C   0  |   | v |   |  d  |
  # which may fail when C is not full rank.

endfunction

