# Quadratic programming
#   min x'*A*x/2 + x'*b
#   s.t. C*x = d
function [x, y] = qprogrameq(A, b, C, d)

  # Reduce the constraints
  #   C*x = d
  # by the transform
  #   x = Z*z + q
  # where C*Z = 0
  [Z, q, Y] = null_space(C, d);

  # The problem is reduced to
  #   min x'*M*x/2 + x'*g
  # where M = Z'*A*Z  and
  #       g = Z'*(b + A*q)
  M = Z'*A*Z;
  g = Z'*(b + A*q);

  # Solve the reduced problem
  z = M\-g;
  x = Z*z + q;

  # Get the Lagrange multipliers
  if nargout > 1
    y = -(C*Y)'\(Y'*(A*x + b));
  endif

  # Note that when C is not full rank, C*Y is a non-square
  # whose inversion becomes solving a least-square problem.

  # It gives the same result as solved
  # by an argumented Lagrangian system
  # | A  C' | * | x | = | -b |
  # | C  0  |   | y |   |  d |
  # which may fail when C is not full rank.
  # [A, C'; C, zeros(rows(C))]\[-b; d]

endfunction

