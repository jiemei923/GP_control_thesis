clc; clear;
close all;
rng(0);
sim_timer = tic;
x0 = [0.5; -0.3; 0.8; -0.6];
trigger_count = 0;

%% Parameters
Tsim = 70;
dt = 0.01;
N = floor(Tsim/dt) + 1;
q_step = 0.001;

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
PB_norm  = norm(P * B, 2);

epsilon_i = sqrt(n) * q_step / 2;
varepsilon = sqrt(nx) * q_step / 2;
Cq = norm(K1,2) * epsilon_i + norm(K2,2) * epsilon_i;

theta = 0.3;
eta_underline = 0.5;
Lf = 2.2;

x_dim = 4;
y_dim = 2;
MaxDataQuantity = 100;

SigmaN = sqrt(1e-6);
SigmaF = 1;
SigmaL = [0.2;0.2;0.2;0.2];

gp = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, ...
                         SigmaN, SigmaF, SigmaL);

%% Initialize
x_set = zeros(nx, N);
x_set(:,1) = x0;

e_tilde_norm_set = zeros(1, N);
e_true_norm_set = zeros(1, N);
z_true_norm_set = zeros(1, N);

trigger_index = zeros(1, N);
trigger_time  = zeros(1, N);
trigger_z_norm = zeros(1, N);

rho_s_set = nan(1, N);
lhs_set = nan(1, N);
rhs_set = nan(1, N);

%% For GP learning performance plot
f_true_set = zeros(y_dim, N);
fhat_gp_set = zeros(y_dim, N);
fhat_used_set = zeros(y_dim, N);
gp_error_norm_set = zeros(1, N);
eta_eval_set = nan(1, N);

for k = 1:N-1

    t_now = (k-1)*dt;
    x_now = x_set(:,k);
    xq_now = quantize_state(x_now, q_step);

    [z_true_now, ~, ~, ~, ~, ~] = compute_z(x_now, t_now);
    z_true_norm_set(k) = norm(z_true_now, 2);

    [zq_now, z1q_now, z2q_now, xd, xd_dot, xd_ddot] = compute_z(xq_now, t_now);

    %% Event-trigger
    do_update_GP = false;

    if k == 1
        do_update_GP = true;
    else
        e_tilde = xk_q - xq_now;
        e_true = xk - x_now;

        e_tilde_norm_set(k) = norm(e_tilde, 2);
        e_true_norm_set(k) = norm(e_true, 2);

        lhs_new = 2 * PB_norm * ...
            (eta_hat_k ...
            + Lf * sqrt(norm(e_tilde, 2) + 2 * varepsilon) ...
            + Lf * sqrt(varepsilon) ...
            + Cq);

        rhs_new = theta * lambda_Qmin * norm(zq_now, 2) + eta_underline;

        rho_s_new = lhs_new - rhs_new;

        lhs_set(k) = lhs_new;
        rhs_set(k) = rhs_new;
        rho_s_set(k) = rho_s_new;

        if rho_s_new >= 0
            do_update_GP = true;
        end
    end

    %% GP update
    if do_update_GP

        trigger_count = trigger_count + 1;
        trigger_index(trigger_count) = k;
        trigger_time(trigger_count) = t_now;
        trigger_z_norm(trigger_count) = z_true_norm_set(k);

        xk = x_now;
        xk_q = quantize_state(xk, q_step);

        fk = f_true(xk);
        yk = fk + SigmaN * randn(size(fk));
        yk_q = quantize_state(yk, q_step);

        if gp.check_Saturation()
            gp.downdateParam(1);
        end

        gp.addPoint(xk_q, yk_q);

        [mu_k, ~, eta_hat_k, ~, ~] = gp.predict(xk_q);
        fhat_k = mu_k;
    end

    %% GP prediction performance, only for plotting
    f_now = f_true(x_now);
    [mu_eval, ~, eta_eval, ~, ~] = gp.predict(xq_now);

    f_true_set(:,k) = f_now;
    fhat_gp_set(:,k) = mu_eval;
    fhat_used_set(:,k) = fhat_k;
    gp_error_norm_set(k) = norm(f_now - mu_eval, 2);
    eta_eval_set(k) = eta_eval;

    %% Controller
    u_now = control_law_gp(z1q_now, z2q_now, xd_ddot, fhat_k, K1, K2);

    %% Simulation
    x_new = system_step(x_now, u_now, dt);
    x_set(:,k+1) = x_new;

end

%% Trim trigger arrays
trigger_index = trigger_index(1:trigger_count);
trigger_time = trigger_time(1:trigger_count);
trigger_z_norm = trigger_z_norm(1:trigger_count);

t_vec = (0:N-1) * dt;

%% Tracking error figure
figure;
semilogy(t_vec, z_true_norm_set);
hold on;
semilogy(trigger_time, trigger_z_norm, 'o', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.2, ...
    'LineStyle', 'none');
xlabel('Time (s)');
ylabel('||z||');
grid on;

figure;

subplot(2,1,1);
plot(t_vec, f_true_set(1,:), 'k-', 'LineWidth', 1.2);
hold on;
plot(t_vec, fhat_used_set(1,:), 'r-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('f_1');
legend('True value', 'GP prediction');
grid on;

subplot(2,1,2);
plot(t_vec, f_true_set(2,:), 'k-', 'LineWidth', 1.2);
hold on;
plot(t_vec, fhat_used_set(2,:), 'r-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('f_2');
legend('True value', 'GP prediction');
grid on;

%% Save GP learning status for comparison
save('result_event_triggered_gpstatus.mat', ...
     't_vec', ...
     'f_true_set', ...
     'fhat_gp_set', ...
     'fhat_used_set', ...
     'gp_error_norm_set', ...
     'eta_eval_set');