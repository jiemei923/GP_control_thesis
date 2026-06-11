clc; clear;
% close all;

%% Go to the folder of this script
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if ~isempty(script_dir)
    cd(script_dir);
end

%% Load results
et = load('result_event_triggered_learning.mat');
tt = load('result_time_triggered_learning.mat');

%% Extract update times
et_time = et.trigger_time(:);
tt_time = tt.update_time(:);

% Remove possible empty values
et_time = et_time(~isnan(et_time));
tt_time = tt_time(~isnan(tt_time));

%% Simulation end time
T_end = et.t_vec(end);

%% Build cumulative update curves
[t_et_step, n_et_step] = build_cumulative_step(et_time, T_end);
[t_tt_step, n_tt_step] = build_cumulative_step(tt_time, T_end);

%% Plot
figure;

plot(t_et_step, n_et_step, 'LineWidth', 1.3);
hold on;
plot(t_tt_step, n_tt_step, 'LineWidth', 1.3);

xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Number of GP updates', 'Interpreter', 'latex');

legend('Event-triggered Learning', ...
       'Time-triggered Learning', ...
       'Location', 'northwest');

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
% This part saves the cumulative update curves into txt files.
% Overleaf will read these txt files and draw the figure.

opt = struct();

% ET data
opt.fname = fullfile(outDir, 'update_count_et_data');
opt.ndata = length(t_et_step);
opt.var_names = {'t','N'};
opt.minval = -1e5;
opt.maxval = 1e5;

data2txt(opt, t_et_step(:), n_et_step(:));

% TT data
opt = struct();
opt.fname = fullfile(outDir, 'update_count_tt_data');
opt.ndata = length(t_tt_step);
opt.var_names = {'t','N'};
opt.minval = -1e5;
opt.maxval = 1e5;

data2txt(opt, t_tt_step(:), n_tt_step(:));

fprintf('Update count data saved to:\n');
fprintf('  %s\n', fullfile(outDir, 'update_count_et_data.txt'));
fprintf('  %s\n', fullfile(outDir, 'update_count_tt_data.txt'));

%% Print update counts
fprintf('Event-triggered GP updates: %d\n', length(et_time));
fprintf('Time-triggered GP updates: %d\n', length(tt_time));

%% Local function
function [t_step, n_step] = build_cumulative_step(update_time, T_end)

    update_time = update_time(:).';
    update_time = update_time(update_time >= 0 & update_time <= T_end);

    num_updates = length(update_time);

    t_step = zeros(1, 2*num_updates + 1);
    n_step = zeros(1, 2*num_updates + 1);

    t_step(1) = 0;
    n_step(1) = 0;

    idx = 2;

    for i = 1:num_updates
        t_step(idx) = update_time(i);
        n_step(idx) = i - 1;
        idx = idx + 1;

        t_step(idx) = update_time(i);
        n_step(idx) = i;
        idx = idx + 1;
    end

    t_step(idx) = T_end;
    n_step(idx) = num_updates;

end