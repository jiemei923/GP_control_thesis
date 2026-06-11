clc; clear;
close all;

%% Load results
no = load('result_no_learning.mat');
off = load('result_offline_learning.mat');
on = load('result_event_triggered_learning.mat');

%% Plot comparison
figure;

semilogy(no.t_vec, no.z_true_norm_set, 'LineWidth', 1.5);
hold on;
semilogy(off.t_vec, off.z_true_norm_set, 'LineWidth', 1.5);
semilogy(on.t_vec, on.z_true_norm_set, 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('||z||');
title('Comparison of No Learning, Offline Learning, and Online Learning');

legend('No Learning', 'Offline Learning', 'Online Learning');

grid on;
hold off;