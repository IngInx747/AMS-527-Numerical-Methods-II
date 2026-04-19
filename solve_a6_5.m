# -u_xx + pi^2*u = 2pi^2*sin(pi*x)
# u(0) = u(1) = 0
function solve_a6_5()

  e0 = 1;
  h0 = 1;

  for i = [10 20 40 80 160 320 640 1280]

    N = i;
    h = 1/N;

    x = linspace(0, 1, N+1);
    x = x(2 : end-1);

    # coefficient matrix
    A = spdiags(repmat([-1, 2 + pi^2*h^2,  -1], N-1, 1), -1:1, N-1, N-1);

    # rhs
    b = sin(pi*x)*pi^2*h^2*2;

    # numerical solution
    u = linsolve(A, b')';

    #plot(x, u); hold on;

    # exact solution
    v = sin(pi*x);

    # error
    e = max(abs(u - v));

    #printf("err = %f\n", e);

    if e0 != 1
      r = log2(e0/e)/log2(h0/h);
      printf("order = %f\n", r);
    endif

    e0 = e;
    h0 = h;

  endfor

endfunction

