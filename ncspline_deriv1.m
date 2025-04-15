# Natural cubic spline: Evaluate 1st derivative
function yy = ncspline_deriv1(xs, s1, xx)

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
        yy(i, j) = s1(1,4) + s1(1,3);
      elseif k >= n
        yy(i, j) = s1(end,3) + s1(end,4);
      else
        a = xval - xs(k);
        b = xs(k+1) - xval;
        h = [a^2; b^2; 1; 1];
        yy(i, j) = s1(k,:)*h;
      endif
    endfor
  endfor

endfunction

