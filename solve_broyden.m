# Broyden's method
function [x, iter] = solve_broyden(f, B, x, tol, max_iter)

  #[x, iter] = solve_broyden_direct(f, B, x, tol, max_iter); return;
  [x, iter] = solve_broyden_smw(f, B, x, tol, max_iter); return;

endfunction

# Broyden's method
function [x, iter] = solve_broyden_direct(f, B, x, tol, max_iter)

  y = f(x);

  for iter = 1 : max_iter
    # update searching direction
    s = B \ -y;
    # move to the new position
    x += s;
    y = f(x);
    if norm(y) < tol && norm(s) < tol
      break
    endif
    # update Jacobian surrogate
    B += y*s' / (s'*s);
  endfor

endfunction

# Broyden's method with Sherman–Morrison–Woodbury formula
function [x, iter] = solve_broyden_smw(f, B, x, tol, max_iter)

  y = f(x);
  H = inv(B);

  for iter = 1 : max_iter
    # update searching direction
    s = -H*y;
    # move to the new position
    x += s;
    y = f(x);
    if norm(y) < tol && norm(s) < tol
      break
    endif
    # update Jacobian surrogate
    s /= s'*s;
    H -= H*y*s'*H / (1 + s'*H*y);
  endfor

endfunction

#printf("(%f, %f, %f, %f)\n", x(1),x(2),x(3),x(4));

