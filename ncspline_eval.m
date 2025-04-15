# Natural cubic spline: Evaluate values
function yy = ncspline_eval(xs, s0, xx)

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
        b = xs(1) - xval;
        y = ncspline_node(xs, s0, 1);
        d = s0(1,3) - s0(1,4);
        yy(i, j) = -d*b + y;
      elseif k >= n
        a = xval - xs(end);
        y = ncspline_node(xs, s0, n);
        d = s0(end,3) - s0(end,4);
        yy(i, j) = d*a + y;
      else
        a = xval - xs(k);
        b = xs(k+1) - xval;
        h = [a^3; b^3; a; b];
        yy(i, j) = s0(k,:)*h;
      endif
    endfor
  endfor

endfunction

# Get y_k on the node_k
function y = ncspline_node(xs, s0, k)

  n = max(size(xs));

  if k >= n # xs = x_n
    d = xs(n) - xs(n - 1);
    h = [d^3; 0; d; 0];
    y = s0(end,:)*h;
  else # xs = x_k
    d = xs(k + 1) - xs(k);
    h = [0; d^3; 0; d];
    y = s0(k,:)*h;
  endif

endfunction

