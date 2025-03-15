# Secant method
function solve_8_1()

  solve_8_1_1();
  solve_8_1_2();

endfunction

# Solve with exact Newtonian method
# f = x^3 - 2x - 5 = 0
function solve_8_1_1()

  f = @(x)(x^3 - x*2 - 5);
  J = @(x)(x^2*3 - 2);
  x = 2.5; # x_0

  tol = 1e-10;
  max_iter = 100000;

  [x, iter] = solve_newton(f, J, x, tol, max_iter);
  y = f(x);

  printf("x_sol = %f\n", x);
  printf("f(x_sol) = %f\n", y);
  printf("iterations: %d\n", iter);

endfunction

# Solve with secant method
# f = x^3 - 2x - 5 = 0
function solve_8_1_2()

  f = @(x)(x^3 - x*2 - 5);
  x0 = 2; x1 = 3;

  tol = 1e-10;
  max_iter = 100000;

  [x, iter] = solve_secant(f, x0, x1, tol, max_iter);
  y = f(x);

  printf("x_sol = %f\n", x);
  printf("f(x_sol) = %f\n", y);
  printf("iterations: %d\n", iter);

endfunction

