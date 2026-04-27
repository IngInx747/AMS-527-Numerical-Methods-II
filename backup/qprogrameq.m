# Quadratic programming
#   min |A*x - b|^2/2
#   s.t. C*x = d
function [x, y] = qprogrameq(A, b, C, d)

  # Reduce the constraints
  #   C*x = d
  # by the transform
  #   x = Z*z + q
  # where C*Z = 0
  [Z, q, Y] = null_space(C, d);

  # The problem is reduced to
  #   min  |M*z - g|^2/2
  # where M = A*Z and
  #       g = b - A*q
  M = A*Z;
  g = b - A*q;

  # Solve the reduced problem
  z = (M'*M)\(M'*g);
  x = Z*z + q;

  # Get the Lagrange multipliers
  if nargout > 1
    y = (C*Y)'\(Y'*A'*(b - A*x));
  endif

  # Note that when C is not full rank, C*Y is a non-square
  # whose inversion becomes solving a least-square problem.

  # It gives the same result as solved
  # by an argumented Lagrangian system
  # | A'A  C' | * | x | = | A'b |
  # |  C   0  |   | y |   |  d  |
  # which may fail when C is not full rank.
  # [A'*A, C'; C, zeros(rows(C))]\[A'*b; d]

endfunction

