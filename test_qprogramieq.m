#
function test_qprogramieq()

  # f = ((x-xc)^2 + (y-yc)^2)/2
  #   = (x^2 + y^2 - 2*xc*x - 2*yc*y + xc^2 + yc^2)/2
  #   = [x y]*I*[x; y]/2 - [x y]*[xc; yc] + b^2/2
  # A = I(2)/2, b = -[xc; yc]
  A = [
    1, 0;
    0, 1];
  b = -[
    1;
    1];
  f = @(x)(x'*A*x + x'*b*2 + b'*b)/2;

  # x <= 1
  # y <= 1
  # x >= 0
  # y >= 0
  # x+y <= 1.5
  # x+y >= 0.5
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
  printf("|f(x_sol)| = %f\n", f(x));
  printf("iterations: %d\n", iter);
  plot(xs(1,:), xs(2,:), ...
    "color", "k", ...
    "marker", ".", ...
    "markersize", 15, ...
    "linestyle", "-", ...
    "linewidth", 2);

endfunction

