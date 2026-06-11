clc; clear;
close all;

%% User settings
full_ndata = 3500;      % 全局图保存的数据点数量，越大越精细，但Overleaf越慢

use_zoom = true;        % true: 生成局部放大图数据；false: 不生成
zoom_t_min = 8;         % 局部放大起始时间
zoom_t_max = 16;        % 局部放大结束时间
zoom_ndata = 1500;      % 局部放大图保存的数据点数量

%% Go to script folder
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if ~isempty(script_dir)
    cd(script_dir);
end

%% Load ET result
et = load('result_event_triggered_learning.mat');

t_vec = et.t_vec(:);
rho_s_set = et.rho_s_set(:);
trigger_time = et.trigger_time(:);

%% Remove initialization update
if numel(trigger_time) >= 2
    trigger_time_plot = trigger_time(2:end);
else
    trigger_time_plot = trigger_time;
end

%% Remove NaN values from trigger function
valid_idx = ~isnan(rho_s_set);
t_phi = t_vec(valid_idx);
rho_s = rho_s_set(valid_idx);

%% Trigger marks are plotted at y = 0
trigger_y = zeros(size(trigger_time_plot));

%% Output folder
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% Full trigger function data
opt = [];
opt.fname = 'figures/trigger_function_data';
opt.ndata = full_ndata;
opt.var_names = {'t','rho_s'};
opt.minval = -1e5;
opt.maxval = 1e5;
data2txt(opt, t_phi, rho_s);

%% Full trigger instant data
opt = [];
opt.fname = 'figures/trigger_time_data';
opt.ndata = length(trigger_time_plot);
opt.var_names = {'t','y'};
opt.minval = -1e5;
opt.maxval = 1e5;
data2txt(opt, trigger_time_plot, trigger_y);

%% Zoomed trigger function data
if use_zoom
    idx_zoom = t_phi >= zoom_t_min & t_phi <= zoom_t_max;
    t_phi_zoom = t_phi(idx_zoom);
    rho_s_zoom = rho_s(idx_zoom);

    idx_trigger_zoom = trigger_time_plot >= zoom_t_min ...
        & trigger_time_plot <= zoom_t_max;
    trigger_time_zoom = trigger_time_plot(idx_trigger_zoom);
    trigger_y_zoom = zeros(size(trigger_time_zoom));

    opt = [];
    opt.fname = 'figures/trigger_function_zoom_data';
    opt.ndata = zoom_ndata;
    opt.var_names = {'t','rho_s'};
    opt.minval = -1e5;
    opt.maxval = 1e5;
    data2txt(opt, t_phi_zoom, rho_s_zoom);

    opt = [];
    opt.fname = 'figures/trigger_time_zoom_data';
    opt.ndata = length(trigger_time_zoom);
    opt.var_names = {'t','y'};
    opt.minval = -1e5;
    opt.maxval = 1e5;
    data2txt(opt, trigger_time_zoom, trigger_y_zoom);
end

%% Print information
fprintf('Generated txt files for TikZ.\n');
fprintf('Full data points: %d\n', full_ndata);
fprintf('Trigger count excluding initial update: %d\n', length(trigger_time_plot));

if use_zoom
    fprintf('Zoom interval: [%.2f, %.2f] s\n', zoom_t_min, zoom_t_max);
    fprintf('Zoom data points: %d\n', zoom_ndata);
    fprintf('Zoom trigger count: %d\n', length(trigger_time_zoom));
end