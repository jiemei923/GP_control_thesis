function [z_now, z1_now, z2_now, xd, xd_dot, xd_ddot] = compute_z(x_now, t_now) %获取期望轨迹，误差

    [xd, xd_dot, xd_ddot] = desired_traj(t_now);

    x1_now = x_now(1:2);
    x2_now = x_now(3:4);

    z1_now = x1_now - xd;
    z2_now = x2_now - xd_dot;
    z_now = [z1_now; z2_now];
end