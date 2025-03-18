# Test Newton method against Rosenbrock function
function test_newton()

  f = @rosenbrock;
  #x0 = [2.5; 2.0];
  x0 = [-1.2; 1.0];

  tol = 1e-6;
  max_iter = 1000;

  hold on;

  [x, iter, xs] = newton(f, x0, tol, max_iter, false);
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

  [x, iter, xs] = newton(f, x0, tol, max_iter, true);
  printf("---- Newton method with line-search ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "r", ...
    "marker", ".", ...
    "markersize", 20, ...
    "linestyle", "--", ...
    "linewidth", 2);

  hold off;

endfunction

