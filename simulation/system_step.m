function x_next = system_step(x, u, dt)

tspan = [0 dt];
[~, x_sol] = ode45(@(t,x) SystemDynamicfun(t, x, u), tspan, x);
x_next = x_sol(end,:)';

end