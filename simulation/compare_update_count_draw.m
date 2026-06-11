clc; clear;
close all;

%% Go to script folder
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if ~isempty(script_dir)
    cd(script_dir);
end

%% User settings
ndata = 1500;

%% Load results
et = load('result_event_triggered_learning.mat');
tt = load('result_time_triggered_learning.mat');

et_time = et.trigger_time(:);
tt_time = tt.update_time(:);

et_time = et_time(et_time >= 0);
tt_time = tt_time(tt_time >= 0);

T_end = min(et.t_vec(end), tt.t_vec(end));

%% Build cumulative update curves
[t_et_step, n_et_step] = build_cumulative_step(et_time, T_end);
[t_tt_step, n_tt_step] = build_cumulative_step(tt_time, T_end);

%% Output folder
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% Save ET cumulative update data
opt = [];
opt.fname = 'figures/update_count_et_data';
opt.ndata = ndata;
opt.var_names = {'t','N'};
opt.minval = -1e5;
opt.maxval = 1e5;
data2txt(opt, t_et_step(:), n_et_step(:));

%% Save TT cumulative update data
opt = [];
opt.fname = 'figures/update_count_tt_data';
opt.ndata = ndata;
opt.var_names = {'t','N'};
opt.minval = -1e5;
opt.maxval = 1e5;
data2txt(opt, t_tt_step(:), n_tt_step(:));

fprintf('Generated cumulative update count txt files.\n');
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