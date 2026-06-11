% function u = control_law_gp(x, fhat_k, K1, K2)
%     x1 = x(1:2);
%     x2 = x(3:4);
%     u = -fhat_k - K1 * x1 - K2 * x2;
% end
function u = control_law_gp(z1, z2, xd_ddot, fhat, K1, K2)
    u = -fhat + xd_ddot - K1*z1 - K2*z2;
end