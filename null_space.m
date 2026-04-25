# Reduce the variables by removing the null space of constraints
#   C*x = d
# The original and the free variables follow the transform
#   x = T*z + q
# where z is the variables of the unconstraint system.
function [T, q, z] = null_space(C, d, x)

  nv = columns(C); # number of variables
  nc = rows(C); # number of constraints

  # C*P = Q*R
  [Q,R,P] = qr(C);

  n0 = rank(R); # number of independent constraints
  n1 = nv - n0; # number of free(reduced) variables

  # Linear independent components can be obtained as
  # C' <- [ R0 | R1 ] * P
  # d' <- (Q' * d)[0 : n0)
  # but they are not needed in the following steps.

  # R = | R0 | R1 |       n0  rows
  #     |  0 |  0 | (nc - n0) rows
  R0 = R(1:n0, 1:n0);
  R1 = R(1:n0, n0+1:end);

  # After decomposition, the constraints become
  # C*x = d => Q*R * P'*x = d
  #         => [ R0 | R1 ]    * P'*x =     Q'*d
  #         => [ I  | R0\R1 ] * P'*x = R0\(Q'*d)
  # In RHS, (Q'*d)[n0+1:nv] should be all zeros.

  # Let v := P'*x = | v0 |  n0 rows
  #                 | v1 |  n1 rows
  # [ I | R0\R1 ] * | v0 | = R0\(Q'*d)
  #                 | v1 |
  # gives v0 = -R0\R1 * v1 + R0\(Q'*d).

  # Hence we have the matrix form
  # | v0 | = | -R0\R1 * v1 + R0\(Q'*d) | =
  # | v1 |   |          v1             |
  # | -R0\R1 | * v1 + | R0\(Q'*d) | = T*v1 + q
  # |    I   |        |     0     |
  # Let z := v1 which is the reduced variables.

  # We can construct q as
  # | R0\(Q'*d) | n0 rows
  # |     0     | n1 rows
  q = Q'*d;
  q = resize(q, n0, 1);
  q = R0\q;
  q = resize(q, nv, 1);

  # We can construct T as
  # | -R0\R1 |  n0 rows
  # |    I   |  n1 rows
  T = [-R0\R1; eye(n1)];

  # Apply permutation from QR
  #   so that x = P*(T*z + q)
  T = P*T;
  q = P*q;

  # Apply transform to the initial values
  # z := v1 := (P'*x)[n0+1:nv]
  if nargout > 2
    z = P'*x;
    z = z(n0+1:end);
  endif

endfunction

