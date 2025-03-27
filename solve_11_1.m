#
function solve_11_1()

  f = @target;
  ce = @equalities;
  ci = @inequalities;
  x0 = [0; 0];

  tol = 1e-10;
  max_iter = 1000;

  printf("---- Quadratic Penalty method ----\n");
  [x, iter, xs] = qpenalty(f, ce, ci, x0, 1., 10., tol, max_iter);
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

function [f, J] = equalities(x)

  f = x(1) + x(2)*2 - 3; # == 0

  if nargout > 1
    J = [1, 2];
  endif

endfunction

function inequalities(x)

  # nothing

endfunction

