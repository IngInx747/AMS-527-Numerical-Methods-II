# Quadratic programming
#   min |A*x - b|^2/2
#   s.t. C*x = d
#        G*x<= h
function [x, iter, xs] = qprogram(A, b, C, d, G, h, x, tol, max_iter)

  # Reduce equality constraints
  #   C*x = d
  # by the transform
  #   x = T*z + q
  [T, q, z] = reducelec(C, d, x);

  # The quadratic optimal problem becomes
  #   min  |A_*z - b_|^2/2
  #   where A_ = A*T
  #         b_ = b - A*q
  # The inequality constraints become
  #   s.t.  G_*z <= h_
  #   where G_ = G*T
  #         h_ = h - G*q
  A_ = A*T;
  b_ = b - A*q;
  G_ = G*T;
  h_ = h - G*q;

endfunction

