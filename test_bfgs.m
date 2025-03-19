# Test BFGS method against Rosenbrock function
function test_bfgs()

  f = @rosenbrock;
  #x0 = [2.5; 2.0];
  x0 = [-1.2; 1.0];

  tol = 1e-6;
  max_iter = 5000;

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
  [x, iter, xs] = bfgsd(f, x0, tol, max_iter);
  printf("---- BFGS method ----\n");
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

  if 0
  [x, iter, xs] = bfgs(f, x0, tol, max_iter);
  printf("---- BFGS method ----\n");
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

  if 1
  [x, iter, xs] = bfgslm(f, x0, 8, tol, max_iter);
  printf("---- Limited-Memory BFGS method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "g", ...
    "marker", ".", ...
    "markersize", 20, ...
    "linestyle", "--", ...
    "linewidth", 2);
  endif

  hold off;

endfunction

