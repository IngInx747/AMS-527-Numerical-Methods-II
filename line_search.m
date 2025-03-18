# Line search algorithm (N.W algo 3.2)
# f: target function and gradient
# x: current position
# s: descent direction dx
# l: maximum step length
# r: reduction factor \in (0, 1) (default 0.5)
# c: Armijo parameter \in (0, 1) (default 1e-4)
# w: Wolfe  parameter \in (0, 1) (default 0.9)
function a = line_search(f, x, s, max_iter, l=20., r=.5, c=1e-4, w=.9)

  [y_0, g_0] = feval(f, x);
  d_0 = g_0'*s;  # initial df
  y_p = y_0 + 1; # last value
  a_p = 0; # last step length
  a = 1;# initial step length

  for iter = 1 : max_iter
    # evaluate function at new step
    [y, g] = feval(f, x + s*a);
    d = g'*s;
    # check Armijo condition
    if y > y_0 + d_0*a*c || y >= y_p
      a = zoomin(f, x, s, y_0, d_0, a_p, a, c, w, max_iter);
      return
    endif
    # check Wolfe condition
    if abs(d) <= -d_0*w
      return
    endif
    # check curvature
    if d >= 0
      a = zoomin(f, x, s, y_0, d_0, a, a_p, c, w, max_iter);
      return
    endif
    # update step length
    a_p = a;
    y_p = y;
    a = a*(1 - r) + l*r;
  endfor

endfunction

# Zoom algorithm (N.W algo 3.3)
# f: target function and gradient
# x: current position
# s: descent direction dx
# [y_0, d_0]: function and df at x = x_0
# [a_1, a_2]: interval of interest
# c: Armijo parameter \in (0, 1) (default 1e-4)
# w: Wolfe  parameter \in (0, 1) (default 0.9)
function a = zoomin(f, x, s, y_0, d_0, a_1, a_2, c, w, max_iter)

  [y_1, g_1] = feval(f, x + s*a_1);

  for iter = 1 : max_iter
    # bisection
    a = (a_1 + a_2)*.5;
    # evaluate function
    [y, g] = feval(f, x + s*a);
    d = g'*s;
    # check Armijo condition
    if y > y_0 + d_0*a*c || y >= y_1
      a_2 = a;
    else
      # check Wolfe condition
      if abs(d) <= -d_0*w
        return
      elseif d*(a_2 - a_1) >= 0
        a_2 = a_1;
      endif
      a_1 = a;
    endif
  endfor

endfunction

