#
function solve_a4_5()

  #verify(); return

  f = @(x)tanh(x*5);
  n = 12; # degree
  n += 1; # points

  x0 = linspace(-1, 1, n);
  y0 = feval(f, x0);
  w0 = bcwgts(x0);

  x1 = cos(pi/2*(2*[0:n-1] + 1)/n);
  y1 = feval(f, x1);
  w1 = bcwgts(x1);

  xs = linspace(-1, 1, 201);

  lerp0 = @(s)bceval(s, x0, y0, w0);
  lerp1 = @(s)bceval(s, x1, y1, w1);

  plot(xs, f(xs), 'color', 'g'); hold on;
  plot(xs, lerp0(xs), 'color', 'r');
  plot(xs, lerp1(xs), 'color', 'b');
  hold off;

endfunction

function verify()

  x = [0 1 2];
  y = [1 4 3];
  w = bcwgts(x);

  x_eval = 0.5;
  y_eval = bceval(x_eval, x, y, w);
  printf("y_eval = %f\n", y_eval);

endfunction

# Barycentric evaluation
function ys = bceval(xs, x, y, w)

  tol = 1e-10;
  ys = zeros(size(xs));

  for k = 1 : max(size(xs))
    v = xs(k) - x;
    I = 0;
    for i = 1 : max(size(x))
      if abs(v(i)) < tol
        I = i;
        break
      endif
    endfor
    if I != 0
      ys(k) = y(i);
    else
      v = w./v;
      ys(k) = sum(v.*y)/sum(v);
    endif
  endfor

endfunction

# Barycentric weights
function w = bcwgts(x)

  n = max (size(x));
  w = ones(size(x));

  for i = 1 : n
    for j = 1 : n
      if i != j
        w(i) *= x(i) - x(j);
      endif
    endfor
  endfor

  w = 1 ./w;

endfunction
