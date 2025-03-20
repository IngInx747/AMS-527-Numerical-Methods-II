# Nonlinear Conjugate Gradient Algorithm
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
function [x, iter, xs] = conjgrad(f, x, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  if nargout(f) < 2
    printf("Target function requires C1 smoothness!\n");
    iter = 0; return
  endif

  # evaluate function
  [~, g] = feval(f, x);
  g_p = g;
  s = - g;

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, 16, :);
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
    # enforce conjugacy
    b_PR = g'*(g - g_p)/(g_p'*g_p); # Polak-Ribiere
    b_HS = g'*(g - g_p)/(s'*(g - g_p)); # Hestenes-Stiefel
    b = max(0, b_HS); # ensure desency
    # searching direction
    s = -g + s*b;
    if norm(s) < tol
      return
    endif
    # backup results
    g_p = g;
  endfor

endfunction

