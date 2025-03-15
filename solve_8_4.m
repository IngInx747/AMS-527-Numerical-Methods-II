# Limited-memory Broyden's method
function solve_8_4()

  K1 = 2;
  K2 = 5;
  A0 = 10;
  B0 = 8;

  f = @(x)[ ...
    x(1)*x(2)*K1 - x(3); ...
    x(3)*K2 - x(4); ...
    x(1) + x(3) + x(4) - A0; ...
    x(2) + x(3) + x(4) - B0];
  J = @(x)[ ...
    x(2)*K1, x(1)*K1, -1, 0; ...
    0, 0, K2, -1; ...
    1, 0, 1, 1; ...
    0, 1, 1, 1];
  x0 = [8; 6; 1; 1];

  if rank(J(x0)) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 1000;
  m = 10;

  #f([7.8551 5.8551 0.9149 2.0851])
  #norm(f([7.8551 5.8551 0.9149 2.0851]))

  [x, iter] = solve_newton(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Newton method ----\n");
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_broyden_limem(f, eye(4), x0, m, tol, max_iter);
  y = f(x);

  printf("---- Limited-memory Broyden's method with B_0 = I ----\n");
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_broyden_limem(f, J(x0), x0, m, tol, max_iter);
  y = f(x);

  printf("---- Limited-memory Broyden's method with B_0 = J(x_0) ----\n");
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

endfunction

