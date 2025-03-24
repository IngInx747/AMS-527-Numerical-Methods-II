#
function solve_11_2()

  f = @target;
  c_eq = @constraint_eq;
  c_ieq = @constraint_ieq;
  x0 = [0; 0];

  tol = 1e-10;
  max_iter = 1000;

  printf("---- Quadratic Penalty method ----\n");
  [x, iter, xs] = quadratic_penalty(f, c_eq, c_ieq, x0, 1., 10., tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f  (x_sol)| = %f\n", norm(f(x)));
  if nargout(c_eq) > 0
    printf("|eq (x_sol)| = %f\n", norm(c_eq(x)));
  endif
  if nargout(c_ieq) > 0
    printf("|ieq(x_sol)| = %f\n", norm(max(c_ieq(x), 0)));
  endif
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "k", ...
    "marker", ".", ...
    "markersize", 15, ...
    "linestyle", "-", ...
    "linewidth", 2);

endfunction

function [f, g, H] = target(x)

  a = 2;
  b = 1;

  f = (x(1) - a)^2 + (x(2) - b)^2;

  if nargout > 1
    g = [ ...
      (x(1) - a)*2; ...
      (x(2) - b)*2];
  endif

  if nargout > 2
    H = [ ...
      2, 0; ...
      0, 2];
  endif

endfunction

function constraint_eq(x)

  # return nothing

endfunction

function [f, J] = constraint_ieq(x)

  f = x(1)^2 + x(2)^2 - 2; # <= 0

  if nargout > 1
    J = [x(1)*2, x(2)*2];
  endif

endfunction

