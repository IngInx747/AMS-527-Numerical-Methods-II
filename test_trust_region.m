# Test Trust-Region method against Rosenbrock function
function test_trust_region()

  f = @rosenbrock;
  #x0 = [2.5; 2.0];
  x0 = [-1.2; 1.0];
  [~, g_0, H_0] = feval(f, x0);

  tol = 1e-6;
  max_iter = 1000;

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

  R0 = norm(H_0'*g_0)*.5;
  [x, iter, xs] = trust_region(f, x0, R0, .25, .75, .5, 2., tol, max_iter);
  printf("---- Trust-region method ----\n");
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

