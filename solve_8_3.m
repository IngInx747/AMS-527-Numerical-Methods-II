# Broyden's method
function solve_8_3()

  f = @the_problem;
  x0 = [8; 6; 1; 1];

  tol = 1e-6;
  max_iter = 1000;

  printf("---- Newton method ----\n");
  [x, iter] = newton(f, x0, tol, max_iter);
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

  printf("---- Broyden method ----\n");
  [x, iter] = broyden(f, x0, tol, max_iter);
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

  printf("---- Limited-memory Broyden method ----\n");
  [x, iter] = broydenlm(f, x0, 10, tol, max_iter);
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

endfunction

function [phi, f, J] = the_problem(x)

  K1 = 2;
  K2 = 5;
  A0 = 10;
  B0 = 8;

  # equations: f = 0
  f = [ ...
    x(1)*x(2)*K1 - x(3); ...
    x(3)*K2 - x(4); ...
    x(1) + x(3) + x(4) - A0; ...
    x(2) + x(3) + x(4) - B0];

  # merit function (required for line-search)
  # Notice: Jacobian is not SPD, line-search problem is not well-defined.
  phi = norm(f);

  if nargout > 2
    # Jacobian
    J = [ ...
      x(2)*K1, x(1)*K1, -1, 0; ...
      0, 0, K2, -1; ...
      1, 0, 1, 1; ...
      0, 1, 1, 1];
  endif

endfunction

function [phi, g, H] = the_problem_opt(x)

  K1 = 2;
  K2 = 5;
  A0 = 10;
  B0 = 8;

  # equations: f = 0
  f = [ ...
    x(1)*x(2)*K1 - x(3); ...
    x(3)*K2 - x(4); ...
    x(1) + x(3) + x(4) - A0; ...
    x(2) + x(3) + x(4) - B0];

  # Jacobian
  J = [ ...
    x(2)*K1, x(1)*K1, -1, 0; ...
    0, 0, K2, -1; ...
    1, 0, 1, 1; ...
    0, 1, 1, 1];

  # merit function (required for line-search)
  phi = f'*f*.5;

  if nargout > 1
    g = J'*f; # gradient of merit function
  endif

  if nargout > 2
    H = J'*J; # Hessian of merit function (approximated)
  endif

endfunction

