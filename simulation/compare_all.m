


clc; clear;
close all;

%% Load results
no = load('result_no_learning.mat');
off = load('result_offline_learning.mat');
on = load('result_event_triggered_learning.mat');
time_result = load('result_time_triggered_learning.mat');

%% Downsample for TikZ export
% Larger skip means fewer points and faster Overleaf compilation.
skip = 10;

idx_no = 1:skip:length(no.t_vec);
idx_off = 1:skip:length(off.t_vec);
idx_on = 1:skip:length(on.t_vec);
idx_time = 1:skip:length(time_result.t_vec);

%% Plot comparison
figure;

semilogy(no.t_vec(idx_no), no.z_true_norm_set(idx_no), 'LineWidth', 1.1);
hold on;
semilogy(off.t_vec(idx_off), off.z_true_norm_set(idx_off), 'LineWidth', 1.1);
semilogy(on.t_vec(idx_on), on.z_true_norm_set(idx_on), 'LineWidth', 1.1);
semilogy(time_result.t_vec(idx_time), time_result.z_true_norm_set(idx_time), 'LineWidth', 1.1);

xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('$\|z\|$', 'Interpreter', 'latex');

legend('No Learning', ...
       'Offline Learning', ...
       'Event-triggered Learning', ...
       'Time-triggered Learning', ...
       'Location', 'best');

grid on;
box on;
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 8);

%% Create output folder
outDir = fullfile(pwd, 'figuretxt');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end


%% Export data for TikZ/PGFPlots

opt = struct();
opt.fname = fullfile(outDir, 'tracking_comparison_data'); %文件名
opt.ndata = 7000;  %数据点数量
opt.var_names = {'t', 'NoGP', 'OffGP', 'ET', 'TT'};
opt.minval = 1e-8;
opt.maxval = 1e8;

data2txt(opt, ...
    no.t_vec(:), ...
    no.z_true_norm_set(:), ...
    off.z_true_norm_set(:), ...
    on.z_true_norm_set(:), ...
    time_result.z_true_norm_set(:));

fprintf('Tracking comparison data saved to: %s\n', ...
    fullfile(outDir, 'tracking_comparison_data.txt'));


% %% Export TikZ only
% outTikz = fullfile(outDir, 'tracking_comparison.tex');
% 
% matlab2tikz(outTikz, ...
%     'width', '\columnwidth', ...
%     'height', '0.65\columnwidth', ...
%     'standalone', false, ...
%     'parseStrings', false, ...
%     'showInfo', false);
% 
% fprintf('TikZ figure saved to: %s\n', outTikz);