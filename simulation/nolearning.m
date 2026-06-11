clc; clear;
close all;
rng(0);
x0 = [0.5; -0.3; 0.8; -0.6];

%%Parameters
Tsim = 70; %仿真时长
dt = 0.01; %仿真间隔
N = floor(Tsim/dt) + 1; %一共有多少个x
q_step = 0.001; %量化精度
nx = 4; %x的总维度

K1 = diag([4, 3]); %用于构造矩阵A
K2 = diag([3, 2.5]); %用于构造矩阵A

y_dim = 2; %dimension of f()

%%initialize condition
x_set = zeros(nx, N); %存储x
x_set(:,1) = x0; %2*rand(4,1) - 1; %x_set(:,1) = [2; -1.5; 0.5; -2.0]; %x初始值
z_true_norm_set = zeros(1, N); %用于存储和期望轨迹的误差z的真实值，画图要用
fhat_k = zeros(y_dim, 1); %fhat一直是0

for k = 1:N-1 %开始循环，对于每一个时刻
    t_now = (k-1)*dt;  %当前时间
    x_now = x_set(:,k); %当前的x真实值
    xq_now = quantize_state(x_now, q_step); %当前x做量化
    [z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now); %把x真实值代入计算，得到z的真实值
    z_true_norm_set(k) = norm(z_true_now, 2); %把真实z的值放进数组
    [zq_now, z1q_now, z2q_now, xd, xd_dot, xd_ddot] = compute_z(xq_now, t_now); %把量化后的x代入计算，得到z_tilde和期望轨迹    
    %% Controller
    % u_now = control_law_gp(x_now, fhat_k, K1, K2); %计算控制率的结果u
    u_now = control_law_gp(z1q_now, z2q_now, xd_ddot, fhat_k, K1, K2);
    %% Simulation
    x_new = system_step(x_now, u_now, dt); %得到下一步的x
    x_set(:,k+1) = x_new; %将下一步的x加入数据集
end
%%draw z
t_vec = (0:N-1) * dt;
figure;
semilogy(t_vec, z_true_norm_set);
title(sprintf('No Learning'));
xlabel('Time (s)');
ylabel('||z||');
grid on;


save('result_no_learning.mat', 't_vec', 'z_true_norm_set');


