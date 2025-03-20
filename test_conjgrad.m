# Test nonlinear CG method against Rosenbrock function
function test_conjgrad()

  f = @rosenbrock;
  #x0 = [1.0; 2.0];
  x0 = [-1.2; 1.0];

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

  printf("---- CG method ----\n");
  [x, iter, xs] = conjgrad(f, x0, tol, max_iter);
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

