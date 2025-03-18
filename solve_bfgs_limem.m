# Limited-Memory BFGS method
function [x, iter] = solve_bfgs_limem(G, B, x, m, tol, max_iter, track=@(~)0)

  if !issymmetric(B)
    printf("BFGS: Input Hessian is not symmetric\n");
    iter = 0; return
  endif

  # initialize memory
  N = rows(x);
  S = zeros(N, 0);
  Y = zeros(N, 0);

  g = G(x); # gradient

  for iter = 1 : max_iter
    # record current position
    track(x);
  endfor

  # record last position
  track(x);

endfunction

function U = enqueue(U, u, m)

  if columns(U) < m
    U = [U, u];
  else
    U(:, 1:m-1) = U(:, 2:m);
    U(:, m) = u;
  endif

endfunction

