#
function solve_a4_7()

  f  = @(x)cos(2*pi*x);
  f1 = @(x)sin(2*pi*x)*-pi*2;
  f2 = @(x)cos(2*pi*x)*-(pi*2)^2;

  xs = linspace(0, 1, 6);
  #xs = chebyshev(0, 1, 15);
  ys = f(xs);

  [s0, s1, s2] = ncspline_init(xs, ys);

  #xx = linspace(-0.2, 1.2, 100);
  xx = linspace(0.1, 0.9, 100);

  yy = ncspline_eval(xs, s0, xx);
  subplot (2, 3, 1);
  plot(xx, yy, xx, f(xx));
  axis('square');
  subplot (2, 3, 4);
  plot(xx, abs(yy - f(xx)));
  axis('square');

  yy = ncspline_deriv1(xs, s1, xx);
  subplot (2, 3, 2);
  plot(xx, yy, xx, f1(xx));
  axis('square');
  subplot (2, 3, 5);
  plot(xx, abs(yy - f1(xx)));
  axis('square');

  yy = ncspline_deriv2(xs, s2, xx);
  subplot (2, 3, 3);
  plot(xx, yy, xx, f2(xx));
  axis('square');
  subplot (2, 3, 6);
  plot(xx, abs(yy - f2(xx)));
  axis('square');

endfunction

function x = chebyshev(a, b, n)

  # (-1, 1)
  x = cos(pi/2*(2*[n-1:-1:0] + 1)/n);

  # (a, b)
  x = (x + 1)/2 * (b - a) + a;

endfunction

