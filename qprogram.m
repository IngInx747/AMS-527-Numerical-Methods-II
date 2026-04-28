# Quadratic programming
#   min x'*A*x/2 + x'*b
#   s.t. C_e*x  = d_e
#        C_i*x <= d_i
function x = qprogram(A, b, Ce, de, Ci, di, x, tol, max_iter)

  # Reduce the constraints
  #   C*x = d
  # by the transform
  #   x = Z*z + q
  # where C*Z = 0
  [Z, q, z] = null_space(Ce, de, x);

  # The problem is reduced to
  #   min x'*A_*x/2 + x'*b_
  # where A_ = Z'*A*Z  and
  #       b_ = Z'*(b + A*q)
  b = Z'*(b + A*q);
  A = Z'*A*Z;

  # The inequality constraints become
  #   s.t.  C_*z <= d_
  #   where C_ = C*Z and
  #         d_ = d - C*q
  di -= Ci*q;
  Ci  = Ci*Z;

  z = qprogramieq(A, b, Ci, di, z, tol, max_iter);
  x = Z*z + q;

endfunction

