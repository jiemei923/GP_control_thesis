clc; clear;
close all;
rng(0);
x0 = [0.5; -0.3; 0.8; -0.6];

%% Parameters
Tsim = 70;
dt = 0.01;
N = floor(Tsim/dt) + 1;
q_step = 0.001;

n  = 2;
nx = 4;

K1 = diag([4, 3]);
K2 = diag([3, 2.5]);

x_dim = 4;
y_dim = 2;
MaxDataQuantity = 400;

SigmaN = sqrt(1e-6);
SigmaF = 1;
SigmaL = [0.2;0.2;0.2;0.2];

gp = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, ...
                         SigmaN, SigmaF, SigmaL);

%% Offline training data generation
% 
% x_min = -1.1 * ones(x_dim, 1);
% x_max =  1.1 * ones(x_dim, 1);
% 
% for i = 1:MaxDataQuantity
% 
%     x_train = x_min + (x_max - x_min) .* rand(x_dim, 1);
%     x_train_q = quantize_state(x_train, q_step);
% 
%     y_train = f_true(x_train) + SigmaN * randn(y_dim, 1);
%     y_train_q = quantize_state(y_train, q_step);
% 
%     gp.addPoint(x_train_q, y_train_q);
% 
% end
N_offline = 400;          % number of pre-collected offline samples

x_min = -1 * ones(x_dim, 1);
x_max =  1 * ones(x_dim, 1);

for i = 1:N_offline

    x_train = x_min + (x_max - x_min) .* rand(x_dim, 1);
    x_train_q = quantize_state(x_train, q_step);

    y_train = f_true(x_train) + SigmaN * randn(y_dim, 1);
    y_train_q = quantize_state(y_train, q_step);

    gp.addPoint(x_train_q, y_train_q);

end

%% Initialize condition
x_set = zeros(nx, N);
x_set(:,1) = x0; %2*rand(4,1) - 1;

z_true_norm_set = zeros(1, N);

%% Simulation loop
for k = 1:N-1

    t_now = (k-1)*dt;
    x_now = x_set(:,k);
    xq_now = quantize_state(x_now, q_step);

    [z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now);
    z_true_norm_set(k) = norm(z_true_now, 2);

    [zq_now, z1q_now, z2q_now, xd, xd_dot, xd_ddot] = compute_z(xq_now, t_now);

    %% GP prediction using fixed offline dataset
    [mu_k, ~, ~, ~, ~] = gp.predict(xq_now);
    fhat_k = mu_k;

    %% Controller
    u_now = control_law_gp(z1q_now, z2q_now, xd_ddot, fhat_k, K1, K2);

    %% Simulation
    x_new = system_step(x_now, u_now, dt);
    x_set(:,k+1) = x_new;

end

%% Draw
t_vec = (0:N-1) * dt;
figure;
semilogy(t_vec, z_true_norm_set);
title(sprintf('Offline Learning'));
xlabel('Time (s)');
ylabel('||z||');
grid on;


save('result_offline_learning.mat', 't_vec', 'z_true_norm_set');