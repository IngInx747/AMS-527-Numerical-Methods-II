#
classdef ncspline

  properties
    s0 = [];
    s1 = [];
    s2 = [];
    xs = [];
  endproperties

  methods

    function s = ncspline(x, y)
      [s0, s1, s2] = ncsplinecoef(x, y);
      xs = x;
    endfunction

    function yy = eval(xx)
      yy = ncsplineval(xs, s0, xx);
    endfunction

    function yy = deriv1(xx)
      ;
    endfunction

    function yy = deriv2(xx)
      ;
    endfunction

  endmethods

endclassdef

# Natural cubic spline setup
function [s0, s1, s2] = ncsplinecoef(x, y)

  if isrow(x)
    x = x';
    y = y';
  endif

  # number of points
  # n + 1

  # number of intervals
  n = rows(x) - 1;

  # intervals
  dx = x(2:n+1) - x(1:n);
  dy = y(2:n+1) - y(1:n);

  h0 = dx(1:n-1);
  h1 = dx(2:n);
  g0 = dy(1:n-1);
  g1 = dy(2:n);

  # linear system (n-1 equations)
  M = [h0/6, (h0 + h1)/3, h1/6];
  A = spdiags(M, -1:1, n-1,n-1);
  b = g1./h1 - g0./h0;

  # bending moments
  m = linsolve(A, b);
  m = [0; m; 0];

  # 0-order coefficients
  s0 = [ ...
    m(2:n+1)./(dx*6), ... # (x - x0)^3
    m(1:n)  ./(dx*6), ... # (x1 - x)^3
    y(2:n+1)./dx - m(2:n+1).*dx/6, ... # (x - x0)
    y(1:n)  ./dx - m(1:n)  .*dx/6];    # (x1 - x)

  # 1-order coefficients
  if nargout > 1
    s1 = [ ...
     -m(2:n+1)./(dx*2), ... # (x - x0)^2
      m(1:n)  ./(dx*2), ... # (x1 - x)^2
      y(2:n+1)./dx - m(2:n+1).*dx/6, ... # 1
     -y(1:n)  ./dx - m(1:n)  .*dx/6];    # 1
  endif

  # 2-order coefficients
  if nargout > 2
    s2 = [ ...
      m(2:n+1)./dx, ... # (x - x0)
      m(1:n)  ./dx];    # (x1 - x)
  endif

endfunction

#
function yy = ncsplineval(x, s0, xx)

  n = max(size(x));
  yy = zeros(size(xx));

  for i = 1 : rows(xx)
    for j = 1 : columns(xx)
      xval = xx(i, j);
      # find the interval
      [~, k] = find(x <= xval);
      if isempty(k)
        k = 0;
      else
        k = k(end);
      endif
      # evaluate within the interval
      if k < 1
        b = x(1) - xval;
        y = ncsplinenode(x, s0, 1);
        d = s0(1,4) - s0(1,3);
        yy(i, j) = d*b + y;
      elseif k >= n
        a = xval - x(end);
        y = ncsplinenode(x, s0, n);
        d = s0(end,3) - s0(end,4);
        yy(i, j) = d*a + y;
      else
        a = xval - x(k);
        b = x(k+1) - xval;
        h = [a^3; b^3; a; b];
        yy(i, j) = s0(k,:)*h;
      endif
    endfor
  endfor

endfunction

#
function y = ncsplinenode(x, s0, k)

  n = max(size(x));

  if k >= n # x = x_n
    d = x(n) - x(n - 1);
    h = [d ^3; 0; d; 0];
    y = s0(end,:)*h;
  else # x = x_k
    d = x(k + 1) - x(k);
    h = [0; d ^3; 0; d];
    y = s0(k,:)*h;
  endif

endfunction

