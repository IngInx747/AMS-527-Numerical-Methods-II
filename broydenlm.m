# Limited-Memory Broyden's Quasi-Newton method
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
function [x, iter, xs] = broydenlm(f, x, m, tol, max_iter, do_line_search = true)

  xs = []; # searching history

  [~, g, B] = feval(f, x); # initial gradient and Hessian

  if rank(B) < rows(x)
    printf("Hessian is degenerated at initial guess!\n");
    iter = 0;
    return
  endif

  H = inv(B);
  H_0 = eye(size(B)); #H;

  # initialize memory
  N = rows(x);
  U = zeros(N, 0);
  V = zeros(N, 0);

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # update searching direction
    s = -H*g;
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
    # update memory
    U = enqueue(U, g, m);
    V = enqueue(V, s/(s'*s), m);
    # update Hessian surrogate
    n = columns(U);
    J = eye(n) + V'*H_0*U;
    if rank(J) < n
      printf("Memory is ill-formed!\n");
      return
    endif
    H = H_0 - (H_0*U)*(J\V')*H_0;
    if rank(H) < rows(x)
      printf("Hessian is degenerated!\n");
      return
    endif
  endfor

endfunction

function U = enqueue(U, u, m)

  if columns(U) < m
    U = [U, u];
  else
    U = [U(:, 2:m), u];
  endif

endfunction

