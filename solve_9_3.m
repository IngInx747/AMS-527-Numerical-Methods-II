# 1D linear search
function solve_9_3()

  f = @(x)(x^3 - x*2 + 2);
  J = @(x)(x^2*3 - 2);
  #x0 = 0;
  #x0 = 1;
  x0 = .5;

  tol = 1e-6;
  max_iter = 100;

  [x, iter] = solve_newton(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Newton method ----\n");
  printf("x_sol = %f\n", x);
  printf("f(x_sol) = %e\n", y);
  printf("iterations: %d\n", iter);

  [x, iter] = solve_newton_backtrack(f, J, x0, tol, max_iter);
  y = f(x);

  printf("---- Backtracking Newton method ----\n");
  printf("x_sol = %f\n", x);
  printf("f(x_sol) = %e\n", y);
  printf("iterations: %d\n", iter);

endfunction
