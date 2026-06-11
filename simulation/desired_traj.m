% function [xd, xd_dot, xd_ddot] = desired_traj(t)
%     xd = [sin(t); cos(t)];
%     xd_dot = [cos(t); -sin(t)];
%     xd_ddot = [-sin(t); -cos(t)];
% end
function [xd, xd_dot, xd_ddot] = desired_traj(t)

    w1 = 0.4;
    w2 = 0.5 * sqrt(2);
    w3 = 0.3;
    w4 = 0.4 * sqrt(3);

    xd = [
        0.7 * sin(w1*t) + 0.2 * sin(w2*t);
        0.6 * cos(w3*t) + 0.2 * sin(w4*t)
    ];

    xd_dot = [
        0.7*w1 * cos(w1*t) + 0.2*w2 * cos(w2*t);
       -0.6*w3 * sin(w3*t) + 0.2*w4 * cos(w4*t)
    ];

    xd_ddot = [
       -0.7*w1^2 * sin(w1*t) - 0.2*w2^2 * sin(w2*t);
       -0.6*w3^2 * cos(w3*t) - 0.2*w4^2 * sin(w4*t)
    ];

end