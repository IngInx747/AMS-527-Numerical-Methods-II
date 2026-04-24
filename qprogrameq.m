# Quadratic programming
#   min |A*x - b|^2/2
#   s.t. C*x = d
function x = qprogrameq(A, b, C, d)

  x = zeros(size(b));

  # Reduce equality constraints
  #   C*x = d
  # by the transform
  #   x = T*z + q
  [T, q] = reducelec(C, d, x);

  # The quadratic optimal problem becomes
  #   min  |A_*z - b_|^2/2
  #   where A_ = A*T
  #         b_ = b - A*q
  b -= A*q;
  A = A*T;

  # Solve the reduced problem
  z = (A'*A) \ (A'*b);
  x = T*z + q;

  # The equations give the same result
  # | A'A  C' | * | x | = | A'b |
  # |  C   0  |   | _ |   |  d  |
  # when C is full rank.

endfunction

