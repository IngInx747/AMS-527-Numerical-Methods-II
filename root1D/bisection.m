# Bisection method: 1D root search
function [x, err] = bisection(f, a, b, tol, max_iter)

  [a, b, err] = expand_interval(f, a, b, abs(b - a)*.5, 2., 1e3);

  if err != 0
    x = 0;
    return
  endif

  for iter = 1 : max_iter
    if b - a < tol
      return
    endif
    x = (a + b)*.5; # mid-point
    if abs(f(x)) < tol
      return
    endif
    if sign(f(a)) != sign(f(x))
      b = x; # root in [a, x]
    else
      a = x; # root in [x, b]
    endif
  endfor

  err = 1;

endfunction

function [a, b, err] = expand_interval(f, a, b, h, s, h_max)

  if a > b
    c = a;
    a = b;
    b = c;
  endif

  c = (a + b)*.5;
  err = 0;

  while sign(f(a))*sign(f(b)) > 0
    if h > h_max
      err = 1;
      return
    endif
    a = c - h;
    b = c + h;
    h *= s;
  endwhile

endfunction
