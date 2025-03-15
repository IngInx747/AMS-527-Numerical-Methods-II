# Broyden's method, limited-memory variants
function [x, iter] = solve_broyden_limem(f, B, x, m, tol, max_iter)

  y = f(x);
  H = Q = inv(B); # Q = H_0 = B_0^{-1}

  # initialize memory
  N = rows(x);
  U = zeros(N, 0);
  V = zeros(N, 0);

  for iter = 1 : max_iter
    # update searching direction
    s = -H*y;
    # move to the new position
    x += s;
    y = f(x);
    if norm(y) < tol && norm(s) < tol
      break
    endif
    # update memory
    U = enqueue(U, y, m);
    V = enqueue(V, s / (s'*s), m);
    # update Jacobian surrogate
    M = columns(U);
    I = eye(M);
    J = I + V'*Q*U;
    if rank(J) < M
      break
    endif
    H = Q - (Q*U)*(J\V')*Q;
  endfor

endfunction

function U = enqueue(U, u, m)

  if columns(U) < m
    U = [U, u];
  else
    U = [U(:, 2:m), u];
  endif

endfunction

