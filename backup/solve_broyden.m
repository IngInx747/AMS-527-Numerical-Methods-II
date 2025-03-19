# Broyden's method with Sherman–Morrison–Woodbury formula
function [x, iter] = solve_broyden(f, B, x, tol, max_iter, track=@(~)0)

  r = f(x); # residue
  H = inv(B); # inversed Jacobian surrogate

  for iter = 1 : max_iter
    # record current position
    track(x);
    # update searching direction
    s = -H*r;
    # move to the new position
    x += s;
    r = f(x);
    if norm(r) < tol && norm(s) < tol
      break
    endif
    # update Jacobian surrogate
    s /= s'*s;
    H -= H*r*s'*H / (1 + s'*H*r);
    if rank(H) < rows(x)
      return
    endif
  endfor

  # record last position
  track(x);

endfunction

# Broyden's method
function [x, iter] = solve_broyden_direct(f, B, x, tol, max_iter, track=@(~)0)

  r = f(x); # residue

  for iter = 1 : max_iter
    # record current position
    track(x);
    # update searching direction
    s = -B\r;
    # move to the new position
    x += s;
    r = f(x);
    if norm(r) < tol && norm(s) < tol
      break
    endif
    # update Jacobian surrogate
    B += r*s' / (s'*s);
    if rank(B) < rows(x)
      return
    endif
  endfor

  # record last position
  track(x);

endfunction

