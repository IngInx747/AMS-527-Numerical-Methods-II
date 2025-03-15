# Solve the nD equations with Newtonian method
function solve_7_1()

  solve_7_1_1();
  solve_7_1_2();

endfunction

# f_1 = x^2 + y^2 - 1 = 0
# f_2 = x - y^3 = 0
function solve_7_1_1()

  f = @(x)[x(1)^2 + x(2)^2 - 1; x(1) - x(2)^3];
  J = @(x)[x(1)*2, x(2)*2; 1, -x(2)^2*3];
  x = [1; 1]; # x_0

  tol = 1e-10;
  max_iter = 100000;

  [x, iter] = solve_newton(f, J, x, tol, max_iter);
  y = f(x);

  printf("(x0, x1) = (%f, %f)\n", x(1), x(2));
  printf("f_1(x0, x1) = %f\n", y(1));
  printf("f_2(x0, x1) = %f\n", y(2));
  printf("iterations: %d\n", iter);

endfunction

# f_1 =  x^2 +  y^2 +  z^2 = 1
# f_2 = 2x^2 +  y^2 - 4z   = 0
# f_3 = 3x^2 - 4y   +  z^2 = 0
function solve_7_1_2()

  f = @(x)[ ...
    x(1)^2   + x(2)^2 + x(3)^2 - 1; ...
    x(1)^2*2 + x(2)^2 - x(3)*4; ...
    x(1)^2*3 - x(2)*4 + x(3)^2];
  J = @(x)[ ...
    x(1)*2, x(2)*2, x(3)*2; ...
    x(1)*4, x(2)*2,     -4; ...
    x(1)*6,     -4, x(3)*2];
  x = [1; 1; 1]; # x_0

  tol = 1e-10;
  max_iter = 100000;

  [x, iter] = solve_newton(f, J, x, tol, max_iter);
  y = f(x);

  printf("(x0, x1, x2) = (%f, %f, %f)\n", x(1), x(2), x(3));
  printf("f_1(x0, x1, x2) = %f\n", y(1));
  printf("f_2(x0, x1, x2) = %f\n", y(2));
  printf("f_3(x0, x1, x2) = %f\n", y(3));
  printf("iterations: %d\n", iter);

endfunction
