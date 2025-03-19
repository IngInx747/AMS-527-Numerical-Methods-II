# f_1 =  x^2 +  y^2 +  z^2 = 1
# f_2 = 2x^2 +  y^2 - 4z   = 0
# f_3 = 3x^2 - 4y   +  z^2 = 0
function solve_7_2()

  f = @the_problem;
  x = [1; 1; 1]; # x_0

  tol = 1e-10;
  max_iter = 1000;

  [x, iter] = newton(f, x, tol, max_iter, false);

  printf("---- Newton method ----\n");
  printf("(x0, x1, x2) = (%f, %f, %f)\n", x(1), x(2), x(3));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

endfunction

function [phi, f, J] = the_problem(x)

  f = [ ...
    x(1)^2   + x(2)^2 + x(3)^2 - 1; ...
    x(1)^2*2 + x(2)^2 - x(3)*4; ...
    x(1)^2*3 - x(2)*4 + x(3)^2];

  # merit function (required for line-search)
  # Notice: Jacobian is not SPD, line-search problem is not well-defined.
  phi = norm(f);

  if nargout > 2
    # Jacobian
    J = [ ...
      x(1)*2, x(2)*2, x(3)*2; ...
      x(1)*4, x(2)*2,     -4; ...
      x(1)*6,     -4, x(3)*2];
  endif

endfunction

