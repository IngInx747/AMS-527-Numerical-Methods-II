# Rosenbrock banana function
function solve_rosenbrock(x, y, a, b)

  a = 1; b = 100;
  f = @(x)((x(1) - a)^2 + (x(2) - x(1)^2)^2*b);
  G = @(x)[(x(1) - a)*2 - (x(2) - x(1)^2)*x(1)*b*4; (x(2) - x(1)^2)*b*2];
  H = @(x)[2 - x(2)*b*4 + x(1)^2*b*12, -x(1)*b*4; -x(1)*b*4, b*2];
  #x0 = [2.5; 2.0];
  x0 = [-1.2; 1.0];

  if rank(H(x0)) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 1000;

  hold on;

  if 1
    tracker([]); # reset tracker
    [x, iter] = solve_newton(G, H, x0, tol, max_iter, @tracker);
    printf("---- Newton method ----\n");
    printf("x_sol = (%f, %f)\n", x(1), x(2));
    printf("|f(x_sol)| = %f\n", norm(f(x)));
    printf("iterations: %d\n", iter);
    plot_path("color", "k", "marker", ".", "markersize", 15, "linestyle", "-.", "linewidth", 3);
  endif

  if 0
    tracker([]); # reset tracker
    [x, iter] = solve_newton_backtrack(G, H, x0, tol, max_iter, @tracker);
    printf("---- Linear-search method ----\n");
    printf("x_sol = (%f, %f)\n", x(1), x(2));
    printf("|f(x_sol)| = %f\n", norm(f(x)));
    printf("iterations: %d\n", iter);
    plot_path("color", "r", "marker", ".", "linestyle", "none");
  endif

  if 0
    tracker([]); # reset tracker
    R0 = norm(H(x0)'*G(x0))*.1;
    [x, iter] = solve_newton_trust_region(G, H, x0, R0, .25, .75, .5, 2., tol, max_iter, @tracker);
    printf("---- Trust-region method ----\n");
    printf("x_sol = (%f, %f)\n", x(1), x(2));
    printf("|f(x_sol)| = %f\n", norm(f(x)));
    printf("iterations: %d\n", iter);
    plot_path("color", "b", "marker", ".", "linestyle", "none");
  endif

  if 0
    tracker([]); # reset tracker
    [x, iter] = solve_broyden(G, H(x0), x0, tol, max_iter, @tracker);
    printf("---- Broyden's method with B_0 = J(x_0) ----\n");
    printf("x_sol = (%f, %f)\n", x(1), x(2));
    printf("|f(x_sol)| = %f\n", norm(f(x)));
    printf("iterations: %d\n", iter);
    plot_path("color", "g", "marker", ".", "linestyle", "none");
  endif

  if 0
    tracker([]); # reset tracker
    [x, iter] = solve_broyden_limem(G, H(x0), x0, 64, tol, max_iter, @tracker);
    printf("---- Limited-memory Broyden's method with B_0 = J(x_0) ----\n");
    printf("x_sol = (%f, %f)\n", x(1), x(2));
    printf("|f(x_sol)| = %f\n", norm(f(x)));
    printf("iterations: %d\n", iter);
    plot_path("color", "c", "marker", ".", "linestyle", "none");
  endif

  if 1
    tracker([]); # reset tracker
    [x, iter] = solve_bfgs(G, H(x0), x0, tol, max_iter, @tracker);
    printf("---- BFGS method with B_0 = J(x_0) ----\n");
    printf("x_sol = (%f, %f)\n", x(1), x(2));
    printf("|f(x_sol)| = %f\n", norm(f(x)));
    printf("iterations: %d\n", iter);
    plot_path("color", "b", "marker", ".", "markersize", 10, "linestyle", "-");
  endif

  hold off;

endfunction

function plot_path(varargin)

  global __xs;

  if isempty(__xs)
    return
  endif

  plot(__xs(1,:), __xs(2,:), varargin{:});

endfunction

