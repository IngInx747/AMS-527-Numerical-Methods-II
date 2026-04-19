# Heun's method
function solve_a5_5()

  f = @(u, t)u*sin(t);
  F = @(t)exp(1 - cos(t));

  t0 = 0;
  t1 = pi*4;

  r = [];
  e0 = 1;
  d0 = 1;

  for i = [10 20 40 80 160 320 640 1280]

    dt = (t1 - t0)/i;
    T = round((t1 - t0)/dt); # num of time steps

    # reset the time bounds (due to rounding error)
    t1 = T*dt;

    u = solve_heun(f, 1, dt, T);

    v = F(t1);
    e = abs(u - v);
    printf("err = %f\n", e);

    if e0 != 1
      r = [r, log2(e0/e)/log2(d0/dt)];
    endif

    e0 = e;
    d0 = dt;

  endfor

  r
  #plot(r);

endfunction

function u = solve_heun(f, u, dt, T)

  for i = 1 : T
    t = (i - 1)*dt;
    # predictor: forward Euler
    v = u + f(u, t)*dt;
    # corrector
    u += (f(u, t) + f(v, t + dt))*dt*.5;
  endfor

endfunction

