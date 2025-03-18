# Inversed Broyden's Quasi-Newton method
# f: scalar, gradient, ~
# B: initial Hessian surrogate
# x: initial guessing position
function [x, iter] = broyden(f, B, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  [~, g] = feval(f, x); # initial gradient
  H = inv(B); # inversed Jacobian surrogate

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # searching direction
    s = -H*g;
    if norm(s) < tol
      break
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, max_iter);
    elseif
      a = 1.;
    endif
    # update position
    x += s*a;
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      break
    endif
    # update Hessian surrogate
    s /= s'*s;
    H -= H*g*s'*H / (1 + s'*H*g);
    if rank(H) < rows(x)
      return
    endif
  endfor

endfunction

# Broyden's Quasi-Newton method
# f: scalar, gradient, ~
# B: initial Hessian surrogate
# x: initial guessing position
function [x, iter] = broyden_direct(f, B, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  [~, g] = feval(f, x); # initial gradient

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # searching direction
    s = -B\g;
    if norm(s) < tol
      break
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, max_iter);
    elseif
      a = 1.;
    endif
    # update position
    x += s*a;
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      break
    endif
    # update Hessian surrogate
    B += g*s' / (s'*s);
    if rank(B) < rows(x)
      return
    endif
  endfor

endfunction

