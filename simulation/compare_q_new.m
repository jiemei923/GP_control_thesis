clc; clear;
close all;

%% Common settings
x0_common = [0.5; -0.3; 0.8; -0.6];

Tsim = 70;
dt = 0.01;
N = floor(Tsim/dt) + 1;

q_step_set = [0.0001, 0.0005, 0.001, 0.003, 0.005];

n  = 2;
nx = 4;

K1 = diag([4, 3]);
K2 = diag([3, 2.5]);

A = [zeros(2), eye(2);
    -K1,      -K2];

B = [zeros(2);
     eye(2)];

Q = eye(nx);
P = lyap(A', Q);

lambda_Qmin = min(eig(Q));
PB_norm = norm(P * B, 2);

lambda_Pmin = min(eig(P));
lambda_Pmax = max(eig(P));

K_norm = norm([K1 K2], 2);

theta = 0.3;

x_dim = 4;
y_dim = 2;

MaxDataQuantity = 100;

SigmaN = sqrt(1e-6);
SigmaF = 1;
SigmaL = [0.2; 0.2; 0.2; 0.2];

%% Storage
z_true_all = zeros(length(q_step_set), N);
trigger_count_all = zeros(length(q_step_set), 1);

bar_wq_all = zeros(length(q_step_set), 1);
bar_wf_all = zeros(length(q_step_set), 1);
bar_w_all = zeros(length(q_step_set), 1);
eta_underline_all = zeros(length(q_step_set), 1);
e_bar_all = zeros(length(q_step_set), 1);

%% Loop over quantization resolutions
for iq = 1:length(q_step_set)

    rng(0);

    q_step = q_step_set(iq);
    trigger_count = 0;

    %% Quantization bounds
    epsilon_i = sqrt(n) * q_step / 2;
    varepsilon = sqrt(nx) * q_step / 2;

    Cq = norm(K1,2) * epsilon_i + norm(K2,2) * epsilon_i;

    %% Theory-consistent parameters from Theorem 1
    bar_wq = varepsilon;

    % Output-label quantization is included because yk_q is used below.
    bar_wy = sqrt(y_dim) * q_step / 2;

    % Aggregated observation/numerical error used in simulation.
    bar_wn = SigmaN + bar_wy;

    % Assumed RKHS norm bounds for simulation.
    bar_f_i = ones(y_dim, 1);

    % Lipschitz bound of the SE kernel used in the proof.
    % The code uses a 4-dimensional ARD SE kernel.
    Lkappa = norm(SigmaF^2 * exp(-0.5) ./ SigmaL);

    % bar_wf from Lemma 1, adapted to the ARD kernel implementation.
    bar_wf = sqrt(2 * bar_wq * sum(bar_f_i.^2) * Lkappa);

    % Aggregated data error bound.
    bar_w = bar_wq + bar_wn + bar_wf;

    % beta from the deterministic GP error bound used in the paper.
    beta_th = sqrt(max(bar_f_i.^2) + MaxDataQuantity);

    % eta underline from Theorem 1.
    eta_underline = sqrt(y_dim) * beta_th * bar_w + bar_wf;

    % L_f = bar_wf / sqrt(bar_wq) in Theorem 1.
    Lf = bar_wf / sqrt(bar_wq);

    xi = theta * lambda_Qmin / (2 * PB_norm);

    % Ultimate tracking error bound e_bar from Theorem 1.
    e_bar = ...
        2 * PB_norm / ((1 - theta) * lambda_Qmin) ...
        * sqrt(lambda_Pmax / lambda_Pmin) ...
        * (bar_wf + eta_underline + K_norm * bar_wq);

    bar_wq_all(iq) = bar_wq;
    bar_wf_all(iq) = bar_wf;
    bar_w_all(iq) = bar_w;
    eta_underline_all(iq) = eta_underline;
    e_bar_all(iq) = e_bar;

    %% GP model
    gp = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, ...
                             bar_w, SigmaF, SigmaL);

    %% Initialize
    x_set = zeros(nx, N);
    x_set(:,1) = x0_common;

    e_tilde_norm_set = zeros(1, N);
    e_true_norm_set = zeros(1, N);
    z_true_norm_set = zeros(1, N);

    for k = 1:N-1

        t_now = (k-1) * dt;
        x_now = x_set(:,k);
        xq_now = quantize_state(x_now, q_step);

        [z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now);
        z_true_norm_set(k) = norm(z_true_now, 2);

        [zq_now, z1q_now, z2q_now, ~, ~, xd_ddot] = compute_z(xq_now, t_now);

        %% Event-trigger
        do_update_GP = false;

        if k == 1
            do_update_GP = true;
        else
            % z_q(t) = x_q(t) - x_{q,k(t)}
            e_tilde = xq_now - xk_q;

            e_true = x_now - xk;

            e_tilde_norm_set(k) = norm(e_tilde, 2);
            e_true_norm_set(k) = norm(e_true, 2);

            % phi(t) from Theorem 1.
            phi_new = ...
                Lf * sqrt(norm(e_tilde, 2) + bar_wq) ...
                + eta_hat_k ...
                - xi * max(norm(zq_now, 2) - bar_wq, 0) ...
                - bar_wf ...
                - eta_underline;

            if phi_new >= 0
                do_update_GP = true;
            end
        end

        %% GP update
        if do_update_GP

            trigger_count = trigger_count + 1;

            xk = x_now;
            xk_q = quantize_state(xk, q_step);

            fk = f_true(xk);
            yk = fk + SigmaN * randn(size(fk));
            yk_q = quantize_state(yk, q_step);

            if gp.check_Saturation()
                gp.downdateParam(1);
            end

            gp.addPoint(xk_q, yk_q);

            [mu_k, var_k, ~, ~, ~] = gp.predict(xk_q);

            % sigma_k = sqrt(trace(Sigma_k)), matching the paper.
            var_k = max(var_k, 0);
            sigma_k = sqrt(sum(var_k));

            eta_hat_k = beta_th * sigma_k + bar_wf;

            fhat_k = mu_k;
        end

        %% Controller
        u_now = control_law_gp(z1q_now, z2q_now, xd_ddot, fhat_k, K1, K2);

        %% Simulation
        x_new = system_step(x_now, u_now, dt);
        x_set(:,k+1) = x_new;

    end

    %% Store final tracking error value
    t_now = (N-1) * dt;
    x_now = x_set(:,N);

    [z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now);
    z_true_norm_set(N) = norm(z_true_now, 2);

    z_true_all(iq,:) = z_true_norm_set;
    trigger_count_all(iq) = trigger_count;

end

%% Plot tracking comparison
t_vec = (0:N-1) * dt;

figure;

for iq = 1:length(q_step_set)
    semilogy(t_vec, z_true_all(iq,:), 'LineWidth', 1.1);
    hold on;

end

xlabel('Time (s)');
ylabel('||z||');

legend(arrayfun(@(q) sprintf('q = %.4g', q), q_step_set, ...
    'UniformOutput', false), ...
    'Location', 'best');

grid on;
box on;

%% Display summary table
result_table = table( ...
    q_step_set(:), ...
    bar_wq_all, ...
    bar_wf_all, ...
    bar_w_all, ...
    eta_underline_all, ...
    e_bar_all, ...
    trigger_count_all, ...
    'VariableNames', { ...
        'q_step', ...
        'bar_wq', ...
        'bar_wf', ...
        'bar_w', ...
        'eta_underline', ...
        'e_bar', ...
        'trigger_count'});

disp(result_table);

%% Export quantization comparison data for TikZ/PGFPlots
outDir = fullfile(pwd, 'figuretxt');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

opt = struct();

opt.fname = fullfile(outDir, 'q_comparison_data_new');

% Export all data points.
% In Overleaf, use "each nth point" to control plotting density.
idx = 1:30:length(t_vec);   % 每隔...个点导出一次
opt.ndata = length(idx);

var_names = cell(1, length(q_step_set) + 1);
var_names{1} = 't';

for iq = 1:length(q_step_set)
    var_names{iq+1} = sprintf('q%d', iq);
end

opt.var_names = var_names;

% Since the plot is semilogy, values should be positive.
opt.minval = 1e-8;
opt.maxval = 1e8;

data_cell = cell(1, length(q_step_set));

for iq = 1:length(q_step_set)
    data_cell{iq} = z_true_all(iq,idx)';
end

data2txt(opt, t_vec(idx)', data_cell{:});

fprintf('Quantization comparison data saved to: %s\n', ...
    fullfile(outDir, 'q_comparison_data_new.txt'));

%% Save MATLAB result
save('result_compare_q.mat', ...
     't_vec', ...
     'q_step_set', ...
     'z_true_all', ...
     'trigger_count_all', ...
     'bar_wq_all', ...
     'bar_wf_all', ...
     'bar_w_all', ...
     'eta_underline_all', ...
     'e_bar_all');