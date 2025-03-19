# Test Trust-region method
function solve_9_5()

  f = @the_problem_opt;
  x0 = [2.0; 2.0];
  [~, g_0, H_0] = feval(f, x0);

  if rank(H_0) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 100;

  hold on;

  [x, iter, xs] = newton(f, x0, tol, max_iter);
  printf("---- Newton method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "k", ...
    "marker", ".", ...
    "markersize", 15, ...
    "linestyle", "-", ...
    "linewidth", 2);

  if 0
  [x, iter, xs] = broyden(f, x0, tol, max_iter);
  printf("---- Broyden method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "r", ...
    "marker", ".", ...
    "markersize", 20, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  R0 = .5; #norm(H_0'*g_0)*.5;
  [x, iter, xs] = trust_region(f, x0, R0, .25, .75, .5, 2., tol, max_iter);
  printf("---- Trust-region method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "b", ...
    "marker", ".", ...
    "markersize", 20, ...
    "linestyle", "--", ...
    "linewidth", 2);

  hold off;

endfunction

function [phi, f, J] = the_problem(x)

  # equations: f = 0
  f = [ ...
    x(1)^2    + x(2) - 2; ...
    exp(x(1)) + x(2) - 1];

  # merit function (required for line-search)
  # Notice: Jacobian is not SPD, line-search problem is not well-defined.
  phi = norm(f);

  if nargout > 2
    # Jacobian
    J = [ ...
      x(1)*2,    1; ...
      exp(x(1)), 1];
  endif

endfunction

function [phi, g, H] = the_problem_opt(x)

  # equations: f = 0
  f = [ ...
    x(1)^2    + x(2) - 2; ...
    exp(x(1)) + x(2) - 1];

  # Jacobian
  J = [ ...
    x(1)*2,    1; ...
    exp(x(1)), 1];

  # merit function (required for line-search)
  phi = f'*f*.5;

  if nargout > 1
    g = J'*f; # gradient of merit function
  endif

  if nargout > 2
    H = J'*J; # Hessian of merit function (approximated)
  endif

endfunction

