# Generic trust region method
# f: [scalar, gradient,  Hessian] (optimal problem)
#    [merit, equations, Jacobian] (solving equations)
# x: initial guess
# R: initial trust region radius
function [x, iter, xs] = trust_region(f, x, R, t1, t2, r1, r2, tol, max_iter)

  xs = []; # searching history

  for iter = 1 : max_iter
    # record current position
    if nargout > 2
      xs = [xs, x];
    endif
    # evaluate function
    [y, g, B] = feval(f, x);
    if norm(g) < tol
      return
    endif
    if rank(B) < rows(x)
      error("Hessian is degenerated!");
    endif
    # find the position that minimizes the model in the region
    #p = dogleg_path(g, B, R);
    p = steihaug_toint(g, B, R, tol, max_iter);
    # construct the approximation model at current position
    m = g'*p + p'*B*p*.5;
    # check if the approximated and the actual reductions meet
    y = feval(f, x + p) - y;
    r = y / m;
    if r < t1 # bad approximation, reject step, shrink region
      R *= r1;
    else # good approximation, accept the step
      x += p;
      if r > t2 && norm(p) >= R # hits the boundary, expand region
        R *= r2;
      endif
      if norm(p) < tol
        return
      endif
    endif
  endfor

endfunction

# Dogleg method to solve trust-region subproblem:
#   min  g'*s + s'*B*s / 2
#   s.t. |s| < R
# g: gradient of merit function
# B: approximated Hessian of merit function
# R: current trust region radius
function s = dogleg_path(g, B, R)

  # Newton step
  SPD = true;
  try L = chol(B);
    p_B = -L\(L'\g);
  catch
    SPD = false;
  end_try_catch

  if SPD && norm(p_B) < R
    s = p_B;
    return
  endif

  # Cauchy point
  p_U = -g * (g'*g)/(g'*B*g);

  if norm(p_U) > R
    s = p_U/norm(p_U)*R;
    return
  endif

  # Dogleg path
  s = boundary(p_U, p_B - p_U, R);

endfunction

# Steihaug-Toint truncated conjugate-gradient method
# to solve trust-region subproblem:
#   min  g'*s + s'*B*s / 2
#   s.t. |s| < R
# g: gradient of merit function
# B: approximated Hessian of merit function
# R: current trust region radius
function s = steihaug_toint(g, B, R, tol, max_iter)

  f = @(x)(g'*x + x'*B*x*.5);
  s = zeros(size(g));
  u = f(s);
  r =  g;
  d = -r;

  for iter = 1 : max_iter
    k = d'*B*d;
    # Check negative curvature
    if k <= 0
      s = boundary(s, d, R);
      return
    endif
    # update position
    a = r'*r /k;
    p = s + d*a;
    # Check radius violation
    if norm(p) >= R
      s = boundary(s, d, R);
      return
    endif
    # Check if model is improved
    v = f(p);
    if u > v # model is improved
      s = p;
    else # model is not improved
      return
    endif
    u = v;
    # update searching direction
    rsq = r'*r;
    r += B*d*a;
    if norm(r) < tol
      s = p;
      return
    endif
    d = d*(r'*r/rsq) - r;
  endfor

endfunction

# Find the position p = a + t*b on the
#   boundary such that |a + t*b| = R.
function p = boundary(a, b, R)

  a2 = a'*a;
  ab = a'*b;
  b2 = b'*b;
  t = (sqrt((R^2 - a2)*b2 + ab^2) - ab)/b2;
  p = a + b*t;

endfunction

