# Limited-Memory BFGS method
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
function [x, iter, xs] = bfgslm(f, x, m, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  if nargout(f) < 2
    printf("Target function requires C1 smoothness!\n");
    iter = 0; return
  endif

  # initialize memory
  N = rows(x);
  S = zeros(N, 0);
  Y = zeros(N, 0);

  # 0-th iteration
  for iter = 1 : 1
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      return
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, -g, 10, :);
    else
      a = 1.;
    endif
    # backup results
    x_p = x;
    g_p = g;
    # update position
    x -= g*a;
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      return
    endif
  endfor

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # secant equation
    s = x - x_p;
    y = g - g_p;
    # update memory
    S = enqueue(S, s, m);
    Y = enqueue(Y, y, m);
    # searching direction
    s = -inv_lm(S, Y, g);
    if norm(s) < tol
      return
    endif
    # determine step length
    if do_line_search
      a = line_search(f, x, s, :);
    else
      a = 1.;
    endif
    # backup results
    x_p = x;
    g_p = g;
    # update position
    x += s*a;
    # evaluate function
    [~, g] = feval(f, x);
    if norm(g) < tol
      return
    endif
  endfor

endfunction

# L-BFGS 2-loop recursion
function q = inv_lm(S, Y, q)

  m = columns(S);
  a = zeros(m,1);

  for i = m : -1 : 1
    a(i) = (S(:,i)'*q)/(S(:,i)'*Y(:,i));
    q -= Y(:,i)*a(i);
  endfor

  q *= S(:,m)'*Y(:,m)/(Y(:,m)'*Y(:,m));

  for i = m : -1 : 1
    b = (Y(:,i)'*q)/(S(:,i)'*Y(:,i));
    q += S(:,i)*(a(i) - b);
  endfor

endfunction

function U = enqueue(U, u, m)

  if columns(U) < m
    U = [U, u];
  else
    U = [U(:, 2:m), u];
  endif

endfunction

