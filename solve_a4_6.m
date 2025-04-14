#
function solve_a4_6()

  hold on;
  test_composite_simpson();
  test_romberg();
  hold off;

endfunction

#
function test_composite_simpson()

  f = @(x)(1 ./(1 + x.^2));
  a = 0; b = 2;
  s0 = atan(b) - atan(a);
  x = [];
  y = [];

  for i = 1 : 10
    n = 2^i;
    s = composite_simpson(f, a, b, n);
    e = abs(s - s0);
    x = [x; n];
    y = [y; e];
  endfor

  x = log2(x);
  y = log2(y);
  #y = -y./x;
  plot(x, y);

endfunction

#
function s = composite_simpson(f, a, b, n)

  h = (b - a)/n;
  x = linspace(a, b, n + 1);
  y = f(x);

  w = ones(size(y));
  for i = 2 : n
    if mod(i, 2)
      w(i) = 2;
    else
      w(i) = 4;
    endif
  endfor

  s = sum(y.*w)*(h/3);

endfunction

#
function test_romberg()

  f = @(x)(1 ./(1 + x.^2));
  a = 0; b = 2;
  s0 = atan(b) - atan(a);
  x = [];
  y = [];

  for i = 1 : 10
    n = 2^i;
    s = romberg(f, a, b, i);
    e = abs(s - s0);
    x = [x; n];
    y = [y; e];
  endfor

  x = log2(x);
  y = log2(y);
  #y = -y./x;
  plot(x, y);

endfunction

#
function s = romberg(f, a, b, n)

  v = zeros(n + 1, 1);

  # initial trapezoidal tableau
  for i = 1 : n + 1
    m = 2 ^(i - 1);
    v(i) = composite_trapezoidal(f, a, b, m);
  endfor

  # extrapolation
  for i = 1 : n + 1
    t = 1/(4 ^i - 1);
    w = [0; v(1 : end-1)];
    v += (v - w)*t;
  endfor

  s = v(end);

endfunction

#
function s = composite_trapezoidal(f, a, b, n)

  h = (b - a)/n;
  x = linspace(a, b, n + 1);
  y = f(x);

  w = ones(size(y)) *2;
  w(1) = 1; w(end) = 1;

  s = sum(y.*w)*(h/2);

endfunction

