clc; clear;
close all;

event_result = load('result_event_triggered_learning.mat');
time_result  = load('result_time_triggered_learning.mat');

figure;

semilogy(event_result.t_vec, event_result.z_true_norm_set, 'LineWidth', 1.5);
hold on;
semilogy(time_result.t_vec, time_result.z_true_norm_set, 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('||z||');
title('Event-triggered vs Time-triggered Online Learning');

legend('Event-triggered online learning', ...
       'Time-triggered online learning', ...
       'Location', 'best');

grid on;
hold off;