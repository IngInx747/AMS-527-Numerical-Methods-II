# Quantum mechanics: spherical well
function solve_a2_6()

  solve_a2_6_1();
  solve_a2_6_2();

endfunction

# ground state
# solution: 2.197824, 1.592484
function solve_a2_6_1()

  s = 3.5;
  f = @(x)[ ...
    x(1) + x(2)*tan(x(1)); ...
    x(1)^2*x(2)^2 - s^2];
  J = @(x)[ ...
    1 + x(2)*sec(x(1))^2, tan(x(1)); ...
    x(1)*x(2)^2*2, x(1)^2*x(2)*2];
  x0 = [2.0; 1.0];

  if rank(J(x0)) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 1000;
  m = 10;

  [x, iter] = solve_newton(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Newton method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_newton_backtrack(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Backtracking Newton method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_broyden(f, J(x0), x0, tol, max_iter);
  y = f(x);

  printf("---- Broyden's method with B_0 = J(x_0) ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_broyden_limem(f, J(x0), x0, m, tol, max_iter);
  y = f(x);

  printf("---- Limited-memory Broyden's method with B_0 = J(x_0) ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

endfunction

# 1st excited state
# solution: 3.298047, 1.061234
function solve_a2_6_2()

  s = 3.5;
  f = @(x)[ ...
    1/(x(1)*tan(x(1))) - 1/x(1)^2 - 1/x(2) - 1/x(2)^2; ...
    x(1)^2*x(2)^2 - s^2];
  J = @(x)[ ...
    -1/(x(1)^2*tan(x(1))) - sec(x(1))^2/(x(1)*tan(x(1)^2)) + 2/x(1)^3, ...
    1/x(2)^2 + 2/x(2)^3; ...
    x(1)*x(2)^2*2, x(1)^2*x(2)*2];
  x0 = [3.2; 1.5];

  if rank(J(x0)) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 1000;
  m = 10;

  [x, iter] = solve_newton(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Newton method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_newton_backtrack(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Backtracking Newton method ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_broyden(f, J(x0), x0, tol, max_iter);
  y = f(x);

  printf("---- Broyden's method with B_0 = J(x_0) ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

  [x, iter] = solve_broyden_limem(f, J(x0), x0, m, tol, max_iter);
  y = f(x);

  printf("---- Limited-memory Broyden's method with B_0 = J(x_0) ----\n");
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(y));
  printf("iterations: %d\n", iter);

endfunction

