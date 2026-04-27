#
function test_qprogramieq()

  # (x-xc)^2 + (y-yc)^2
  A = [
    1, 0;
    0, 1];
  b = [
    2;
    2];

  # x <= 1
  # y <= 1
  # x >= 0
  # y >= 0
  # x+y <= 1.5
  Ci = [
    1, 0;
    0, 1;
   -1, 0;
    0,-1;
    1, 1];
  di = [
    1;
    1;
    0;
    0;
    1.5];

  max_iter = 1000;
  tol = 1e-10;
  x = [0; 0];

  printf("---- Quadratic Programming: Active Set method ----\n");
  [x, iter, xs] = qprogramieq(A, b, Ci, di, x, tol, max_iter);

  printf("x_sol = (%f, %f)\n", x(1), x(2));
  printf("|f(x_sol)| = %f\n", norm(A*x-b));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "k", ...
    "marker", ".", ...
    "markersize", 15, ...
    "linestyle", "-", ...
    "linewidth", 2);

endfunction

