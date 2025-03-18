# Broyden–Fletcher–Goldfarb–Shanno method
function [x, iter] = solve_bfgs(G, B, x, tol, max_iter, track=@(~)0)

  if !issymmetric(B)
    printf("BFGS: Input Hessian is not symmetric\n");
    iter = 0; return
  endif

  g = G(x); # gradient
  H = inv(B); # inversed Jacobian surrogate

  for iter = 1 : max_iter
    # record current position
    track(x);
    # update searching direction
    s = -H*g;
    # move to the new position
    x += s;
    y = G(x) - g; # delta gradient
    if norm(y) < tol && norm(s) < tol
      break
    endif
    # update Hessian surrogate
    R = y*s'/(y'*s); # curvature: s'*B*s
    H += R'*H*R - R'*H - H*R + s*s'/(y'*s);
    if rank(H) < rows(x)
      return
    endif
    # update gradient
    g += y;
  endfor

  # record last position
  track(x);

endfunction

# Broyden–Fletcher–Goldfarb–Shanno method
function [x, iter] = solve_bfgs_direct(G, B, x, tol, max_iter, track=@(~)0)

  if !issymmetric(B)
    printf("BFGS: Input Hessian is not symmetric\n");
    iter = 0; return
  endif

  g = G(x); # gradient

  for iter = 1 : max_iter
    # record current position
    track(x);
    # update searching direction
    s = -B\g;
    # move to the new position
    x += s;
    y = G(x) - g; # delta gradient
    if norm(y) < tol && norm(s) < tol
      break
    endif
    # update Hessian surrogate
    B += y*y'/(s'*y) + g*g'/(s'*g);
    if rank(B) < rows(x)
      return
    endif
    # update gradient
    g += y;
  endfor

  # record last position
  track(x);

endfunction

