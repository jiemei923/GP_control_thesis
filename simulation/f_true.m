function f_true = f_true(x)
    x11 = x(1);
    x12 = x(2);
    x21 = x(3);
    x22 = x(4);

    f_true = [1 - sin(x11) + 0.5 / (1 + exp(-x21 / 10));
             0.8 * cos(x12) + 0.3 * x22^3];
    
end