#
function fcontour(f, x, y, varargin)

  [X,Y] = meshgrid(x, y);
  Z = zeros(size(X));

  for i = 1 : rows(X)
    for j = 1 : columns(X)
      Z(i, j) = f([X(i, j); Y(i, j)]);
    endfor
  endfor

  contour(X, Y, Z, varargin{:});

endfunction

