# Muller's method
# Eq: Z(w) = R + iwL + 1/iwC
function solve_8_2()

  f = @(x)(100 + x*.01i - 1e-6i/x);
  x0 = (1 + 1i)*0.5;
  x1 = (1 + 1i)*1.0;
  x2 = (1 + 1i)*1.5;

  tol = 1e-10;
  max_iter = 100000;

  [x, iter] = solve_muller(f, x0, x1, x2, tol, max_iter);
  y = f(x);

  printf("x_sol = %f + %fi\n", real(x), imag(x));
  printf("f(x_sol) = %f + %fi\n", real(y), imag(y));
  printf("iterations: %d\n", iter);

endfunction

