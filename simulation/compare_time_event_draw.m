clc; clear;
close all;

%% Go to script folder
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if ~isempty(script_dir)
    cd(script_dir);
end

%% User settings
ndata = 1500;   % 保存点数，越大越精细，Overleaf越慢

%% Load data
et = load('result_event_triggered_learning.mat');
tt = load('result_time_triggered_learning.mat');

t_et = et.t_vec(:);
z_et = et.z_true_norm_set(:);

t_tt = tt.t_vec(:);
z_tt = tt.z_true_norm_set(:);

%% Remove invalid values
valid_et = ~isnan(z_et) & z_et > 0;
valid_tt = ~isnan(z_tt) & z_tt > 0;

t_et = t_et(valid_et);
z_et = z_et(valid_et);

t_tt = t_tt(valid_tt);
z_tt = z_tt(valid_tt);

%% Output folder
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% Save ET data
opt = [];
opt.fname = 'figures/tracking_et_data';
opt.ndata = ndata;
opt.var_names = {'t','z'};
opt.minval = -1e5;
opt.maxval = 1e5;
data2txt(opt, t_et, z_et);

%% Save TT data
opt = [];
opt.fname = 'figures/tracking_tt_data';
opt.ndata = ndata;
opt.var_names = {'t','z'};
opt.minval = -1e5;
opt.maxval = 1e5;
data2txt(opt, t_tt, z_tt);

fprintf('Generated ET/TT tracking txt files.\n');