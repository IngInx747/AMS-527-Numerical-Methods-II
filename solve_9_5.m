# Linear search v.s. Trust region
function solve_9_5()

  f = @(x)[ ...
    x(1)^2    + x(2) - 2; ...
    exp(x(1)) + x(2) - 1];
  J = @(x)[ ...
    x(1)*2,    1; ...
    exp(x(1)), 1];
  x0 = [2.0; 2.0];

  if rank(J(x0)) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 100;

  [x, iter] = solve_newton_backtrack(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Linear-search method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  R0 = norm(J(x0)'*f(x0))*.1;
  [x, iter] = solve_newton_trust_region(f, J, x0, R0, .25, .75, .5, 2., tol, max_iter);
  y = f(x);

  printf("---- Trust-region method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

endfunction

