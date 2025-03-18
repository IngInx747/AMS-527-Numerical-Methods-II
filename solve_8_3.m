# Broyden's method
function solve_8_3()

  f = @the_problem;
  x0 = [8; 6; 1; 1];
  [~, f_0, J_0] = feval(f, x0);

  if rank(J_0) < rows(x0)
    printf("Jacobian is degenerated at initial guess!\n");
    return
  endif

  tol = 1e-6;
  max_iter = 1000;

  [x, iter] = newton(f, x0, tol, max_iter);
  printf("---- Newton method ----\n");
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

  [x, iter] = broyden(f, eye(4), x0, tol, max_iter, false);
  printf("---- Broyden's method with B_0 = I ----\n");
  printf("x_sol = (%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
  printf("|f(x_sol)| = %f\n", norm(f(x)));
  printf("iterations: %d\n", iter);

  [x, iter] = broyden(f, J_0, x0, tol, max_iter);
  printf("---- Broyden's method with B_0 = J(x_0) ----\n");
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
  phi = f'*f*.5;

  if nargout > 2
    # Jacobian
    J = [ ...
      x(2)*K1, x(1)*K1, -1, 0; ...
      0, 0, K2, -1; ...
      1, 0, 1, 1; ...
      0, 1, 1, 1];
  endif

endfunction

