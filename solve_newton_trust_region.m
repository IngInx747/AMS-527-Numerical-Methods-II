# Generic trust region method
# f: function R^{n} -> R^{n}
# J: Jacobian of f
# x: initial guess
# R: initial trust region radius
function [x, iter] = solve_newton_trust_region(f, J, x, R, t1, t2, r1, r2, tol, max_iter)

  # merit function u(x) = |f(x)|^2/2
  u = @(x)(f(x)'*f(x)*.5);

  # model function m(p) = u(x) + G(x)'*p + p'*H(x)*p/2
  G = @(x)(J(x)'*f(x)); # gradient of u
  H = @(x)(J(x)'*J(x)); # Hessian of u (approximated)

  for iter = 1 : max_iter
    # construct the approximation model
    g = G(x);
    B = H(x);
    if rank(B) < rows(x)
      return
    endif
    # find the path that minimizes the model within the region
    p = dogleg(x, g, B, R);
    # check if the approximated and the actual reductions meet
    y = u(x);
    if norm(sqrt(y*2)) < tol
      return
    endif
    m = y + g'*p + p'*B*p*.5;
    r = (u(x + p) - y)/(m - y);
    if r < t1 # bad approximation, reject step, shrink region
      R *= r1;
    else # good approximation, accept the step
      x += p;
      if r > t2 && norm(p) >= R # step hits the boundary, expand region
        R *= r2;
      endif
      if norm(p) < tol
        return
      endif
    endif
  endfor

endfunction

# x: current position
# g: gradient of merit function
# B: approximated Hessian of merit function
# R: current trust region radius
function p = dogleg(x, g, B, R)

  # Newton step
  SPD = true;
  try L = chol(B);
    p_B = -L\(L'\g);
  catch
    SPD = false;
  end_try_catch

  if SPD && norm(p_B) < R
    p = p_B;
    return
  endif

  # Cauchy point
  p_U = -g * (g'*g)/(g'*B*g);

  if norm(p_U) > R
    p = p_U/norm(p_U)*R;
    return
  endif

  # Dogleg path
  a2 = p_U'*p_U;
  ab = (p_B - p_U)'*p_U;
  b2 = (p_B - p_U)'*(p_B - p_U);
  r = (sqrt((R^2 - a2)*b2 + ab^2) - ab)/b2;
  p = p_U + (p_B - p_U)*r;

endfunction

