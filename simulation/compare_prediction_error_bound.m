clc; clear;
close all;
rng(0);
x0 = 2*rand(4,1) - 1;

%%Parameters
Tsim = 30; %仿真时长
dt = 0.01; %仿真间隔
N = floor(Tsim/dt) + 1; %一共有多少个x
q_step_list = [0.0001, 0.001, 0.01, 0.1, 0.2]; %不同量化精度
eta_time_all = nan(length(q_step_list), N); %存储eta数据

for q_idx = 1:length(q_step_list)
    q_step = q_step_list(q_idx);




n  = 2; %每一阶x的维度
nx = 4; %x的总维度

K1 = diag([4, 3]); %用于构造矩阵A
K2 = diag([3, 2.5]); %用于构造矩阵A

A = [zeros(2), eye(2);
    -K1,      -K2]; %用于生成矩阵P,用于后续触发计算
B = [zeros(2);
     eye(2)];  %用于后续触发计算

Q = eye(nx);
P = lyap(A', Q);

lambda_Qmin = min(eig(Q)); %用于触发计算
PB_norm  = norm(P * B, 2); %用于触发计算
epsilon_i = sqrt(n) * q_step / 2; 
varepsilon = sqrt(nx) * q_step / 2;
Cq = norm(K1,2) * epsilon_i + norm(K2,2) * epsilon_i;

theta = 0.3;
eta_underline = 0.2;
Lf    = 3.0;

x_dim = 4; %x的总维度rn
y_dim = 2; %dimension of f()
MaxDataQuantity = 100; %Dataset最大数量
SigmaN = sqrt(1e-6);
SigmaF = 1; %用于GP
SigmaL = [2;2;2;2]; %用于GP
gp = LocalGP_MultiOutput(x_dim,y_dim,MaxDataQuantity, ...
				SigmaN,SigmaF,SigmaL);


%%initialize condition
x_set = zeros(nx, N); %存储x
x_set(:,1) = x0; %x_set(:,1) = [2; -1.5; 0.5; -2.0]; %x初始值
e_tilde_norm_set = zeros(1, N); %用于存储e_tilde_norm
e_true_norm_set = zeros(1, N); %用于存储e真实值，后面画图要用
z_true_norm_set = zeros(1, N); %用于存储和期望轨迹的误差z的真实值，画图要用

    eta_time_set = nan(1, N);

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
        lhs_new = 2 * PB_norm * (eta_hat_k + Lf * sqrt(norm(e_tilde, 2) + 2 * varepsilon) + Lf * sqrt(varepsilon) + Cq);
        rhs_new = theta * lambda_Qmin * norm(zq_now, 2) + eta_underline;
        if lhs_new >= rhs_new
            do_update_GP = true; %计算不等式的左右两边判断是否满足触发条件
        end
    end
    %% GP
    if do_update_GP %如果触发了
 %first trigger sample
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

    eta_time_set(k) = eta_hat_k;
end

    eta_time_set(N) = eta_time_set(N-1);
    eta_time_all(q_idx, :) = eta_time_set;

end
%%draw画图：横坐标时间-纵坐标log10||z||
t_vec = (0:N-1) * dt;

figure;

for q_idx = 1:length(q_step_list)
    semilogy(t_vec, eta_time_all(q_idx, :), 'LineWidth', 1.2);
    hold on;
end

xlabel('Time (s)');
ylabel('Prediction error bound \eta');

legend_str = strings(1, length(q_step_list));
for q_idx = 1:length(q_step_list)
    legend_str(q_idx) = "q = " + num2str(q_step_list(q_idx));
end
legend(legend_str, 'Location', 'best');

grid on;

% %%draw画图：横坐标时间-纵坐标log10||e||
% t_vec = (0:N-1) * dt;
% figure;
% log10_e_true_norm_set = log10(e_true_norm_set);
% plot(t_vec, log10_e_true_norm_set);
% xlabel('Time (s)');
% ylabel('log10||e||');
% grid on;
% 


