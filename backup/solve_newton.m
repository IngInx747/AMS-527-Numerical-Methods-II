# Newton method - root search
function [x, iter] = solve_newton(f, J, x, tol, max_iter, track=@(~)0)

  # J(x_k)*s_k = -f(x_k)
  # x_{k+1} = x_k + s_k

  for iter = 1 : max_iter
    # record current position
    track(x);
    # evaluate residue
    y = f(x);
    if norm(y) < tol
      break
    endif
    # descent direction
    H = J(x);
    if rank(H) < rows(x)
      return
    endif
    s = H \ -y;
    # update position
    x += s;
    if norm(s) < tol
      break
    endif
  endfor

  # record last position
  track(x);

endfunction
