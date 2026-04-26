# Reduce the variables by keeping the null space of constraints
#   C*x = d
# The original and the free variables(x_z) follow the transform
#   x = Z*x_z + Y*x_y
#     = Z*x_z + q
# where C*Z = 0 and C*Y is full rank.
function varargout = null_space(C, d, varargin)

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
  # | -R0\R1 | * v1 + | R0\(Q'*d) | = Z*v1 + q
  # |    I   |        |     0     |
  # Let x_z := v1 which are the free variables.

  # We can construct Z as
  # | -R0\R1 |  n0 rows
  # |    I   |  n1 rows
  Z = [-R0\R1; eye(n1)];

  # We can construct Y as
  # |  R0\I  |  n0 rows
  # |    0   |  n1 rows
  #Y = resize(R0\eye(n0), nv, n0);

  # We only need C*Y to be full rank so any Y works
  Y = resize(eye(n0), nv, n0);

  # We can construct q as
  # | R0\(Q'*d) | n0 rows
  # |     0     | n1 rows
  x_y = resize(Q'*d, n0, 1);
  q = resize(R0\x_y, nv, 1); # q = Y*x_y

  # Apply permutation to the transform
  #   x = P*(Z*x_z + Y*x_y)
  Z = P*Z;
  Y = P*Y;
  q = P*q;

  # Apply permutation to the initial x
  #   x_z := v1 := (P'*x)[n0+1:nv]
  if nargin > 2
    x_z = (P'*varargin{1})(n0+1:end);
  endif

  # export everything
  varargout{1} = Z;
  varargout{2} = q;

  if nargout == 3
    if nargin > 2 # [Z, q, x_z]
      varargout{3} = x_z;
    else # [Z, q, Y]
      varargout{3} = Y;
    endif
  elseif nargout == 4 # [Z, q, Y, x_z]
    varargout{3} = Y;
    varargout{4} = x_z;
  endif

endfunction

