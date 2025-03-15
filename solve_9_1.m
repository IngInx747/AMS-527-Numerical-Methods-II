# 1D linear search
function solve_9_1()

  f = @(x)(x^3 - x*2 + 2);

  tol = 1e-6;
  max_iter = 100;

  [x, err] = solve_bisection(f, 0, 1, tol, max_iter);
  y = f(x);

  printf("---- Bisection method ----\n");
  printf("x_sol = %f\n", x);
  printf("f(x_sol) = %e\n", y);

endfunction
