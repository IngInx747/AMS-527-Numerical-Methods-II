# Newton method with backtracking linear search
function [x, iter] = solve_newton_backtrack(f, J, x, tol, max_iter, track=@(~)0)

  phi = @(x)(f(x)'*f(x)*.5); # merit function
  grd = @(x)(J(x)'*f(x)); # gradient of f at x
  rho = 0.5;
  c_1 = 1e-4;
  c_2 = .9;
  a_0 = 1.;

  for iter = 1 : max_iter
    # record current position
    track(x);
    # evaluate residue
    y = f(x); # f_k = f(x_k)
    if norm(y) < tol
      return
    endif
    # descent direction
    H = J(x);
    if rank(H) < rows(x)
      return
    endif
    s = H \ -y;
    # determine step length
    [a, err] = backtrack(phi, grd, x, s, a_0, rho, c_1, c_2, max_iter);
    if err != 0
      return
    endif
    # update position
    x += s*a;
    if norm(s) < tol
      break
    endif
  endfor

  # record last position
  track(x);

endfunction

# f: merit function
# g: gradient of f
# x: current position
# s: descent direction
# a: initial step length
# r: reduction factor \in (0, 1) (default 0.5)
# c: Armijo parameter \in (0, 1) (default 1e-4)
# w: Wolfe parameter  \in (0, 1) (default 0.9)
function [a, err] = backtrack(f, g, x, s, a, r, c, w, max_iter)

  f_x = f(x);
  g_x = g(x);
  err = 0;

  for iter = 1 : max_iter
    if f(x + s*a) - f_x <= (g_x'*s)*a*c && ...
       abs(g(x + s*a)'*s) <= abs(g_x'*s)*w
      return
    endif
    a *= r;
  endfor

  err = 1;

endfunction

