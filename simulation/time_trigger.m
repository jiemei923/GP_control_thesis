clc; clear;
close all;
rng(0);
x0 = [0.5; -0.3; 0.8; -0.6];
sim_timer = tic;
trigger_count = 0;%计算触发次数


%%Parameters
Tsim = 70; %仿真时长
dt = 0.01; %仿真间隔
N = floor(Tsim/dt) + 1; %一共有多少个x
q_step = 0.001; %量化精度

n  = 2; %每一阶x的维度
nx = 4; %x的总维度

K1 = diag([4, 3]); %用于构造矩阵A
K2 = diag([3, 2.5]); %用于构造矩阵A

time_update_period = 0.2; % time-trigger更新周期%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_update_steps = round(time_update_period / dt);

x_dim = 4; %x的总维度rn
y_dim = 2; %dimension of f()
MaxDataQuantity = 100; %Dataset最大数量
SigmaN = sqrt(1e-6);
SigmaF = 1; %用于GP
SigmaL = [0.2;0.2;0.2;0.2]; %用于GP
gp = LocalGP_MultiOutput(x_dim,y_dim,MaxDataQuantity, ...
				SigmaN,SigmaF,SigmaL);


%%initialize condition
x_set = zeros(nx, N); %存储x
x_set(:,1) = x0; %2*rand(4,1) - 1; %x_set(:,1) = [2; -1.5; 0.5; -2.0]; %x初始值
e_tilde_norm_set = zeros(1, N); %用于存储e_tilde_norm
e_true_norm_set = zeros(1, N); %用于存储e真实值，后面画图要用
z_true_norm_set = zeros(1, N); %用于存储和期望轨迹的误差z的真实值，画图要用
update_time = zeros(1, N);   % 存储time-triggered GP更新时刻

for k = 1:N-1 %开始循环，对于每一个时刻
    t_now = (k-1)*dt;  %当前时间
    x_now = x_set(:,k); %当前的x真实值
    xq_now = quantize_state(x_now, q_step); %当前x做量化
    [z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now); %把x真实值代入计算，得到z的真实值
    z_true_norm_set(k) = norm(z_true_now, 2); %把真实z的值放进数组
    [zq_now, z1q_now, z2q_now, xd, xd_dot, xd_ddot] = compute_z(xq_now, t_now); %把量化后的x代入计算，得到z_tilde和期望轨迹

    %% Event-trigger
    do_update_GP = false;
    if k == 1
        do_update_GP = true; %判断是否运行do_update_GP,如果第一次则必须运行
    else
        e_tilde = xk_q - xq_now; %计算当前的e_tilde
        e_true = xk - x_now;  %计算当前e的真实值
        e_tilde_norm_set(k) = norm(e_tilde, 2); %把每个时刻e_tilde的值存进e_tilde_norm_set
        e_true_norm_set(k) = norm(e_true,2); %把每个时刻e_true的值放进e_true_norm_set
        if mod(k-1, time_update_steps) == 0
            do_update_GP = true; %计算不等式的左右两边判断是否满足触发条件
        end
    end
    %% GP
    if do_update_GP %如果触发了
 %first trigger sample
        trigger_count = trigger_count + 1; %触发次数+1
        update_time(trigger_count) = t_now;
        xk = x_now; %当前的x
        xk_q = quantize_state(xk, q_step); %获取量化后的x后续用于加入数据集
        fk = f_true(xk);
        yk = fk + SigmaN * randn(size(fk)); %yk = f_true(xk);
        yk_q = quantize_state(yk, q_step); %计算需要放进数据集的y
        if gp.check_Saturation()
            gp.downdateParam(1); %判断是否需要删掉之前数据防止数据集爆满
        end
        gp.addPoint(xk_q,yk_q); %把当前触发点的数据加入数据集
        [mu_k, ~, eta_hat_k, ~, ~] = gp.predict(xk_q);
        fhat_k = mu_k; %计算新的f_hat
    end
    %% Controller
    % u_now = control_law_gp(x_now, fhat_k, K1, K2); %计算控制率的结果u
    u_now = control_law_gp(z1q_now, z2q_now, xd_ddot, fhat_k, K1, K2);
    %% Simulation
    x_new = system_step(x_now, u_now, dt); %得到下一步的x
    x_set(:,k+1) = x_new; %将下一步的x加入数据集
end

    %%draw画图：横坐标时间-纵坐标log10||z||
t_vec = (0:N-1) * dt;
figure;
semilogy(t_vec, z_true_norm_set);
title(sprintf('Time-triggered online learning, update period \\Delta t = %.5f s', time_update_period));
xlabel('Time (s)');
ylabel('||z||');
grid on;

update_time = update_time(1:trigger_count);

save('result_time_triggered_learning.mat', ...
     't_vec', 'z_true_norm_set', ...
     'update_time', 'trigger_count');

fprintf('trigger_count = %.4f ', trigger_count); %输出显示触发次数

elapsed_time = toc(sim_timer);
fprintf('Elapsed time = %.4f s\n', elapsed_time);

