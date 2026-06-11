function dx = SystemDynamicfun(t, x, u)
dx = zeros(4,1);
x1 = x(1:2);
x2 = x(3:4);
dx(1:2) = x2;
dx(3:4) = f_true(x) + u;
end