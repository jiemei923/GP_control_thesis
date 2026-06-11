clc; clear;
close all;
rng(0);
x0_common = [0.5; -0.3; 0.8; -0.6]; % x0_common = 2*rand(4,1) - 1;

%%Parameters
Tsim = 70; %仿真时长
dt = 0.01; %仿真间隔
N = floor(Tsim/dt) + 1; %一共有多少个x
q_step_set = [0.0001, 0.0005, 0.001, 0.003, 0.005,]; %不同量化精度
z_true_all = zeros(length(q_step_set), N); %用于存储所有q对应的z
trigger_count_all = zeros(length(q_step_set), 1); %存储每个q对应的触发次数




for iq = 1:length(q_step_set)
    % rng(0);
    q_step = q_step_set(iq);
    trigger_count = 0;


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

theta = 0.3; %几乎没影响触发次数
eta_underline = 0.5; %越小触发越多
Lf    = 2.2; %越大触发越多

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
x_set(:,1) = x0_common; %x_set(:,1) = [2; -1.5; 0.5; -2.0]; %x初始值
e_tilde_norm_set = zeros(1, N); %用于存储e_tilde_norm
e_true_norm_set = zeros(1, N); %用于存储e真实值，后面画图要用
z_true_norm_set = zeros(1, N); %用于存储和期望轨迹的误差z的真实值，画图要用

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
        trigger_count = trigger_count + 1;
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

t_now = (N-1)*dt;
x_now = x_set(:,N);
[z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now);
z_true_norm_set(N) = norm(z_true_now, 2);

z_true_all(iq,:) = z_true_norm_set;
trigger_count_all(iq) = trigger_count;


end

%%draw画图：横坐标时间-纵坐标log10||z||

t_vec = (0:N-1) * dt;

skip = 1;
idx = 1:skip:length(t_vec);

figure;

% semilogy(t_vec, z_true_all(1,:));
% hold on;
% semilogy(t_vec, z_true_all(2,:));
% semilogy(t_vec, z_true_all(3,:));
% semilogy(t_vec, z_true_all(4,:));
% semilogy(t_vec, z_true_all(5,:));
semilogy(t_vec(idx), z_true_all(1,idx), 'LineWidth', 1.1);
hold on;
semilogy(t_vec(idx), z_true_all(2,idx), 'LineWidth', 1.1);
semilogy(t_vec(idx), z_true_all(3,idx), 'LineWidth', 1.1);
semilogy(t_vec(idx), z_true_all(4,idx), 'LineWidth', 1.1);
semilogy(t_vec(idx), z_true_all(5,idx), 'LineWidth', 1.1);


xlabel('Time (s)');
ylabel('||z||');

legend(arrayfun(@(q) sprintf('q = %.4g', q), q_step_set, ...
    'UniformOutput', false), ...
    'Location', 'best');


grid on;

result_table = table(q_step_set(:), trigger_count_all, ...
    'VariableNames', {'q_step', 'trigger_count'});

disp(result_table);

%% Export quantization comparison data for TikZ/PGFPlots
% This part only saves the data into a txt file.
% Overleaf will read this txt file and draw the figure.

outDir = fullfile(pwd, 'figuretxt');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

opt = struct();

% File name without ".txt"
opt.fname = fullfile(outDir, 'q_comparison_data');

% Export all data points.
% In Overleaf, use "each nth point" to control plotting density.
opt.ndata = length(t_vec);

% Column names for LaTeX.
% The first column is time.
% The following columns correspond to different q_step values.
opt.var_names = {'t', 'q1', 'q2', 'q3', 'q4', 'q5'};

% Since the plot is semilogy, values should be positive.
opt.minval = 1e-8;
opt.maxval = 1e8;


data2txt(opt, ...
    t_vec(:), ...
    z_true_all(1,:)', ...
    z_true_all(2,:)', ...
    z_true_all(3,:)', ...
    z_true_all(4,:)', ...
    z_true_all(5,:)');

fprintf('Quantization comparison data saved to: %s\n', ...
    fullfile(outDir, 'q_comparison_data.txt'));