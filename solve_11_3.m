#
function solve_11_3()

  f = @target;
  ce = @equalities;
  ci = @inequalities;
  x0 = [0; 0];

  tol = 1e-10;
  max_iter = 1000;

  printf("---- Augmented Lagrangian method ----\n");
  [x, iter, xs] = lagrange(f, ce, ci, x0, 100., tol, max_iter);
  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f (x_sol)| = %f\n", norm(f(x)));
  if nargout(ce) > 0
    printf("|Ce(x_sol)| = %f\n", norm(ce(x)));
  endif
  if nargout(ci) > 0
    printf("|Ci(x_sol)| = %f\n", norm(max(ci(x), 0)));
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

  a = 3;
  b = 2;

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

function [f, J] = equalities(x)

  f = x(1) + x(2) - 4; # == 0

  if nargout > 1
      J = [1, 1];
  endif

endfunction

function [f, J] = inequalities(x)

  f = x(1) - 1; # <= 0

  if nargout > 1
    J = [1, 0];
  endif

endfunction

