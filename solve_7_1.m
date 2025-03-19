# f_1 = x^2 + y^2 - 1 = 0
# f_2 = x - y^3 = 0
function solve_7_1()

  f = @the_problem;
  x0 = [1; 1];

  tol = 1e-10;
  max_iter = 1000;

  [x, iter] = newton(f, x0, tol, max_iter, false);

  printf("---- Newton method ----\n");
  printf("(x0, x1) = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

endfunction

function [phi, f, J] = the_problem(x)

  f = [ ...
    x(1)^2 + x(2)^2 - 1; ...
    x(1) - x(2)^3];

  # merit function (required for line-search)
  # Notice: Jacobian is not SPD, line-search problem is not well-defined.
  phi = norm(f);

  if nargout > 2
    # Jacobian
    J = [ ...
      x(1)*2, x(2)*2; ...
      1, -x(2)^2*3];
  endif

endfunction

