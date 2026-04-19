# Runge–Kutta method
function solve_a6_4()

  f = @(u, t)-50*(u - cos(t));
  F = @(t)((50*cos(t) + sin(t))*50 + exp(t*-50))/2501;

  t0 = 0;
  t1 = 1;

  dt = .05;
  T = round((t1 - t0)/dt); # num of time steps

  # reset the time bounds (due to rounding error)
  t1 = T*dt;

  # mesh
  t = linspace(t0, t1, T + 1);

  # RK4
  u1 = solve_RK4(f, 1, dt, T);

  # implicit trapzoidal method
  u2 = solve_IT2(f, 1, dt, T);

  # analytic solution
  v = F(t);

  abs(u1(end) - v(end))
  abs(u2(end) - v(end))

  if 0
  plot(t, u1, 'linewidth', 2, ...
       t, u2, 'linewidth', 2, ...
       t, v, 'k', 'linestyle', '--');
  endif

  if 1
  plot(t, abs(u1 - v), 'linewidth', 2, ...
       t, abs(u2 - v), 'linewidth', 2);
  endif

endfunction

#
function U = solve_RK4(f, u, dt, T)

  U = [u];

  for i = 1 : T
    t = (i - 1)*dt;
    v0 = f(u, t);
    v1 = f(u + v0*dt*.5, t + dt*.5);
    v2 = f(u + v1*dt*.5, t + dt*.5);
    v3 = f(u + v2*dt, t + dt);
    u += (v0 + v1*2 + v2*2 + v3)*dt/6.;
    U = [U, u];
  endfor

endfunction

#
function U = solve_IT2(f, u, dt, T)

  U = [u];

  for i = 1 : T
    t = (i - 1)*dt;
    v0 = f(u, t);
    du = @(x)(x - u - (v0 + f(x, t + dt))*dt*.5);
    u = fsolve(du, u);
    U = [U, u];
  endfor

endfunction

