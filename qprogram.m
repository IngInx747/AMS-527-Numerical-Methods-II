# Quadratic programming
#   min |A*x - b|^2/2
#   s.t. C_e*x  = d_e
#        C_i*x <= d_i
function x = qprogram(A, b, Ce, de, Ci, di, x, tol, max_iter)

  # Reduce the constraints
  #   C*x = d
  # by the transform
  #   x = Z*z + q
  # where C*Z = 0
  [Z, q, z] = null_space(Ce, de, x);

  # The problem becomes
  #   min  |A_*z - b_|^2/2
  #   where A_ = A*Z and
  #         b_ = b - A*q
  b -= A*q;
  A  = A*Z;

  # The inequality constraints become
  #   s.t.  C_*z <= d_
  #   where C_ = C*Z and
  #         d_ = d - C*q
  di -= Ci*q;
  Ci  = Ci*Z;

  z = qprogramieq(A, b, Ci, di, z, tol, max_iter);
  x = Z*z + q;

endfunction

