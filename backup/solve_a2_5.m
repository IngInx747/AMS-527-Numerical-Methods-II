#
function solve_a2_5()

  x0 = 11/2;
  x1 = 61/11;
  x = zeros(1, 0);

  for iter = 1 : 20
    x2 = 111 - (1130 - 3000/x0) / x1;
    x0 = x1;
    x1 = x2;
    x = [x(1,:), x2];
  endfor

  plot(x);

endfunction
