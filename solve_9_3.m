# Example of why linear-search matters in Newton method
function solve_9_3()

  f = @the_problem;

  # pure Newton method stuck in [0, 1]
  x0 = 0;
  #x0 = 1;
  #x0 = .5;

  tol = 1e-6;
  max_iter = 1000;

  [x, iter] = newton(f, x0, tol, max_iter, false);
  printf("---- Newton method ----\n");
  printf("x_sol = %f\n", x);
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

  [x, iter] = newton(f, x0, tol, max_iter);
  printf("---- Newton method with backtrack ----\n");
  printf("x_sol = %f\n", x);
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

endfunction

function [phi, f, J] = the_problem(x)

  # a carefully constructed function
  f = x^3 - x*2 + 2;

  if nargout > 2
    J = x^2*3 - 2;
  endif

  # merit function (required for line-search)
  phi = abs(f);

endfunction

