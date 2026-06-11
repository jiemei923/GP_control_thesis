clc; clear;
close all;

%% Go to the folder of this script
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if ~isempty(script_dir)
    cd(script_dir);
end

%% Load results
et = load('result_event_triggered_learning.mat');
tt = load('result_time_triggered_learning.mat');

%% Extract trigger/update times
et_time = et.trigger_time(:);
tt_time = tt.update_time(:);

et_time = et_time(~isnan(et_time));
tt_time = tt_time(~isnan(tt_time));

%% Remove possible empty entries
et_time = et_time(et_time >= 0);
tt_time = tt_time(tt_time >= 0);

%% Plot trigger instances
figure;

hold on;

plot(et_time, 2 * ones(size(et_time)), 'x', ...
    'Color', [0.0000 0.4470 0.7410], ...
    'MarkerSize', 5, ...
    'LineWidth', 1.0);

plot(tt_time, 1 * ones(size(tt_time)), 'x', ...
    'Color', [0.8500 0.3250 0.0980], ...
    'MarkerSize', 5, ...
    'LineWidth', 1.0);

yticks([1 2]);
yticklabels({'Time-triggered', 'Event-triggered'});

xlabel('Time (s)', 'Interpreter', 'latex');

xlim([0 max(et.t_vec)]);
ylim([0.5 2.5]);

grid on;
box on;

set(gcf, 'Color', 'w');
set(gca, 'FontSize', 8);

%% Create output folder
outDir = fullfile(pwd, 'figures');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Export TikZ
outTikz = fullfile(outDir, 'trigger_instances_comparison.tex');

matlab2tikz(outTikz, ...
    'width', '\columnwidth', ...
    'height', '0.45\columnwidth', ...
    'standalone', false, ...
    'parseStrings', false, ...
    'showInfo', false);

fprintf('TikZ figure saved to: %s\n', outTikz);

%% Print counts
fprintf('Event-triggered updates: %d\n', length(et_time));
fprintf('Time-triggered updates: %d\n', length(tt_time));