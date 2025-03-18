# rosenbrock.m This function returns the function value, partial derivatives
# and Hessian of the (general dimension) rosenbrock function, given by:
#
#   f(x) = sum_{i=1:N-1} 100*(x_{i+1} - x_i^2)^2 + (1 - x_i)^2
#
# where N is the dimension of x. The true minimum is 0 at x = (1 1 ... 1).
#
# Carl Edward Rasmussen, 2001-07-21.

function [f, G, H] = rosenbrock(x)

  a = 1;
  b = 100;

  N = length(x);
  f = sum((x(2:N) - x(1:N-1).^2).^2*b + (a - x(1:N-1)).^2);

  if nargout > 1
    G = zeros(N, 1);
    G(1:N-1) = -x(1:N-1).*(x(2:N) - x(1:N-1).^2)*b*4 - (a - x(1:N-1))*2;
    G(2:N) += (x(2:N) - x(1:N-1).^2)*b*2;
  endif

  if nargout > 2
    H = zeros(N, N);
    H(1:N-1, 1:N-1) = diag(-x(2:N)*b*4 + x(1:N-1).^2*b*12 + 2);
    H(2:N, 2:N) += eye(N-1)*b*2;
    H -= diag(x(1:N-1)*b*4, 1) + diag(x(1:N-1)*b*4, -1);
  endif

endfunction

