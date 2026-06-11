clc; clear;
close all;

%% Go to script folder
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if ~isempty(script_dir)
    cd(script_dir);
end

%% Load ET and TT GP learning status
et = load('result_event_triggered_gpstatus.mat');
tt = load('result_time_triggered_learning.mat');

%% Use common time length
N_common = min(length(et.t_vec), length(tt.t_vec)) - 1;

t = et.t_vec(1:N_common).';

% True value is taken along the ET closed-loop trajectory
f1_true = et.f_true_set(1,1:N_common).';
f2_true = et.f_true_set(2,1:N_common).';

% GP predictions actually used by ET and TT controllers
f1_et = et.fhat_used_set(1,1:N_common).';
f2_et = et.fhat_used_set(2,1:N_common).';

f1_tt = tt.fhat_used_set(1,1:N_common).';
f2_tt = tt.fhat_used_set(2,1:N_common).';

%% Output folder
outDir = fullfile(pwd, 'figures');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Export data
opt = struct();
opt.fname = fullfile(outDir, 'gp_prediction_comparison_data');
opt.ndata = length(t);     % export all data points
opt.var_names = {'t', 'f1_true', 'f1_ET', 'f1_TT', ...
                      'f2_true', 'f2_ET', 'f2_TT'};
opt.minval = -1e5;
opt.maxval = 1e5;

data2txt(opt, ...
    t, ...
    f1_true, f1_et, f1_tt, ...
    f2_true, f2_et, f2_tt);

fprintf('GP prediction comparison data saved to: %s\n', ...
    fullfile(outDir, 'gp_prediction_comparison_data.txt'));