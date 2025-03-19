# Broyden's Quasi-Newton method (direct Hessian)
# f: scalar, gradient, Hessian
# x: initial guessing position
function [x, iter, xs] = broydend(f, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  if nargout(f) >= 3
    [~, g, B] = feval(f, x); # initial gradient and Hessian
  elseif nargout(f) >= 2
    [~, g] = feval(f, x); # initial gradient
    B = eye(rows(x));
  else
    printf("Target function requires C1 smoothness!\n");
    iter = 0; return
  endif

  if rank(B) < rows(x)
    printf("Initial Hessian is degenerated!\n");
    iter = 0; return
  endif

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # searching direction
    s = -B\g;
    if norm(s) < tol
      return
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, :);
    else
      a = 1.;
    endif
    # update position
    x += s*a;
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      return
    endif
    # update Hessian surrogate
    B += g*s' / (s'*s);
    if rank(B) < rows(x)
      printf("Hessian is degenerated!\n");
      return
    endif
  endfor

endfunction

