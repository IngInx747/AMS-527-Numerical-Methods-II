# Quadratic programming
#   min |A*x - b|^2/2
#   s.t. C*x = d
function [x, v] = qprogrameq(A, b, C, d)

  # Reduce equality constraints
  #   C*x = d
  # by the transform
  #   x = Z*z + y
  # where C*Z = 0
  [Z, y, Y] = null_space(C, d);

  # The quadratic optimal problem becomes
  #   min  |M*z - g|^2/2
  #   where M = A*Z
  #         g = b - A*y
  M = A*Z;
  g = b - A*y;

  # Solve the reduced problem
  z = (M'*M)\(M'*g);
  x = Z*z + y;

  # Lagrange multipliers
  if nargout > 1
    v = (C*Y)'\(Y'*A'*(b - A*x));
  endif

  # This method gives the same result as solved
  # by an argumented Lagrangian system
  # | A'A  C' | * | x | = | A'b |
  # |  C   0  |   | v |   |  d  |
  # which may fail when C is not full rank.

endfunction

