# Rosenbrock banana function
function solve_rosenbrock(x, y, a, b)

  a = 1; b = 100;
  f = @(x)((x(1) - a)^2 + (x(2) - x(1)^2)^2*b);
  G = @(x)[(x(1) - a)*2 - (x(2) - x(1)^2)*x(1)*b*4; (x(2) - x(1)^2)*b*2];
  H = @(x)[2 - x(2)*b*4 + x(1)^2*b*12, -x(1)*b*4; -x(1)*b*4, b*2];
  x0 = [2.0; 2.0];

  if rank(H(x0)) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 1000;

  [x, iter] = solve_newton(G, H, x0, tol, max_iter);
  y = f(x);

  printf("---- Newton method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_newton_backtrack(G, H, x0, tol, max_iter);
  y = f(x);

  printf("---- Linear-search method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  R0 = norm(H(x0)'*G(x0))*.1;
  [x, iter] = solve_newton_trust_region(G, H, x0, R0, .25, .75, .5, 2., tol, max_iter);
  y = f(x);

  printf("---- Trust-region method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

endfunction

