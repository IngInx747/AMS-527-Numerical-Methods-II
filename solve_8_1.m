# Quasi-Newton method
# f = x^3 - 2x - 5 = 0
function solve_8_1()

  f = @the_problem;
  x0 = 2.5;
  x1 = 3;

  tol = 1e-10;
  max_iter = 1000;

  printf("---- Newton method ----\n");
  [x, iter] = newton(f, x0, tol, max_iter, false);
  printf("x_sol = %f\n", x);
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

  printf("---- Secant method ----\n");
  [x, iter] = secant(f, x0, x1, tol, max_iter);
  printf("x_sol = %f\n", x);
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

endfunction

function [phi, f, J] = the_problem(x)

  f = x^3 - x*2 - 5;

  if nargout > 2
    J = x^2*3 - 2;
  endif

  # merit function (required for line-search)
  phi = 0; # doesn't matter since we don't use line-search

endfunction

