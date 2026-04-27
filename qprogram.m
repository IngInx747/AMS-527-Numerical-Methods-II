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
  #   min x'*M*x/2 + x'*g
  # where M = Z'*A*Z  and
  #       g = Z'*(b + A*q)
  M = Z'*A*Z;
  g = Z'*(b + A*q);

  # The inequality constraints become
  #   s.t.  C_*z <= d_
  #   where C_ = C*Z and
  #         d_ = d - C*q
  di -= Ci*q;
  Ci  = Ci*Z;

  z = qprogramieq(M, g, Ci, di, z, tol, max_iter);
  x = Z*z + q;

endfunction

