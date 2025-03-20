# Quantum mechanics: spherical well
# ground state (2.197824, 1.592484)
function solve_a2_6_0()

  f = @the_problem_opt;
  x0 = [2.0; 2.0];

  tol = 1e-6;
  max_iter = 1000;

  hold on;

  printf("---- Newton method ----\n");
  [x, iter, xs] = newton(f, x0, tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "k", ...
    "marker", ".", ...
    "markersize", 15, ...
    "linestyle", "-", ...
    "linewidth", 2);

  if 1
  printf("---- Trust-region method ----\n");
  R0 = .5; #norm(H_0'*g_0)*.5;
  [x, iter, xs] = trust_region(f, x0, R0, .25, .75, .5, 2., tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "b", ...
    "marker", ".", ...
    "markersize", 20, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  if 1
  printf("---- Broyden method ----\n");
  [x, iter, xs] = broyden(f, x0, tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "r", ...
    "marker", "o", ...
    "markersize", 5, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  if 1
  printf("---- Limited-memory Broyden method ----\n");
  [x, iter, xs] = broydenlm(f, x0, 8, tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "r", ...
    "marker", "x", ...
    "markersize", 5, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  if 1
  printf("---- BFGS method ----\n");
  [x, iter, xs] = bfgs(f, x0, tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "g", ...
    "marker", "o", ...
    "markersize", 5, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  if 1
  printf("---- Limited-Memory BFGS method ----\n");
  [x, iter, xs] = bfgslm(f, x0, 8, tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "g", ...
    "marker", "x", ...
    "markersize", 5, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  hold off;

endfunction

function [phi, f, J] = the_problem(x)

  s = 3.5;

  # equations: f = 0
  f = [ ...
    x(1) + x(2)*tan(x(1)); ...
    x(1)^2*x(2)^2 - s^2];

  # merit function (required for line-search)
  # Notice: Jacobian is not SPD, line-search problem is not well-defined.
  phi = norm(f);

  if nargout > 2
    # Jacobian
    J = [ ...
      1 + x(2)*sec(x(1))^2, tan(x(1)); ...
      x(1)*x(2)^2*2, x(1)^2*x(2)*2];
  endif

endfunction

function [phi, g, H] = the_problem_opt(x)

  s = 3.5;

  # equations: f = 0
  f = [ ...
    x(1) + x(2)*tan(x(1)); ...
    x(1)^2*x(2)^2 - s^2];

  # Jacobian
  J = [ ...
    1 + x(2)*sec(x(1))^2, tan(x(1)); ...
    x(1)*x(2)^2*2, x(1)^2*x(2)*2];

  # merit function (required for line-search)
  phi = f'*f*.5;

  if nargout > 1
    g = J'*f; # gradient of merit function
  endif

  if nargout > 2
    H = J'*J; # Hessian of merit function (approximated)
  endif

endfunction

