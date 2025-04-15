# Natural cubic spline: Evaluate 2nd derivative
function yy = ncspline_deriv2(xs, s2, xx)

  n = max(size(xs));
  yy = zeros(size(xx));

  for i = 1 : rows(xx)
    for j = 1 : columns(xx)
      xval = xx(i, j);
      # find the interval
      [~, k] = find(xs <= xval);
      if isempty(k)
        k = 0;
      else
        k = k(end);
      endif
      # evaluate within the interval
      if k < 1
        yy(i, j) = 0;
      elseif k >= n
        yy(i, j) = 0;
      else
        a = xval - xs(k);
        b = xs(k+1) - xval;
        h = [a; b];
        yy(i, j) = s2(k,:)*h;
      endif
    endfor
  endfor

endfunction

