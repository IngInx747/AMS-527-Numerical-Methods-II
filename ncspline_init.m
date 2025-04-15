# Natural cubic spline: Setup coefficients table
function [s0, s1, s2] = ncspline_init(x, y)

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
      m(2:n+1)./(dx*2), ... # (x - x0)^2
     -m(1:n)  ./(dx*2), ... # (x1 - x)^2
      y(2:n+1)./dx - m(2:n+1).*dx/6, ... # 1
     -y(1:n)  ./dx + m(1:n)  .*dx/6];    # 1
  endif

  # 2-order coefficients
  if nargout > 2
    s2 = [ ...
      m(2:n+1)./dx, ... # (x - x0)
      m(1:n)  ./dx];    # (x1 - x)
  endif

endfunction

