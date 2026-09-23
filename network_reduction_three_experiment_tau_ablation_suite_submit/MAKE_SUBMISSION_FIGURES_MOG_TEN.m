%% MAKE_SUBMISSION_FIGURES_MOG_TEN.M
%
% Postprocessing only.
%
% PURPOSE
% -------
% This script reads the COMPLETED results of the three experiments:
%
%   results_01_noise_tau_ablation_...
%   results_02_q_tau_ablation_...
%   results_03_sigma_tau_ablation_...
%
% and recreates the three submission figures WITHOUT rerunning any solver.
%
% The new figures:
%   1) use "MoG-TeN" instead of "Proposed";
%   2) use thicker curves / markers / axes;
%   3) use larger, bold axis labels and tick labels;
%   4) keep the same experimental data and numerical results;
%   5) are saved under:
%
%        submission_figures_MoG-TeN/
%
%
% OUTPUT FILES
% ------------
%   submission_noise_sweep_X_NMSE.png
%   submission_q_sweep_X_NMSE.png
%   submission_sigma_sweep_X_NMSE.png
%
% Each figure is also saved as .fig and vector .pdf.
%
%
% EXPECTED LOCATION
% -----------------
% Place this script in the SAME folder as the three RUNME files and the
% completed results folders, for example:
%
%   network_reduction_three_experiment_tau_ablation_suite/
%       RUNME_1_noise_sweep_with_tau_ablation.m
%       RUNME_2_q_sweep_with_tau_ablation.m
%       RUNME_3_sigma_sweep_with_tau_ablation.m
%       results_01_noise_tau_ablation_.../
%       results_02_q_tau_ablation_.../
%       results_03_sigma_tau_ablation_.../
%       MAKE_SUBMISSION_FIGURES_MOG_TEN.m
%


clear; clc; close all;

%% ========================================================================
% 0. LOCATE THIS SCRIPT
% =========================================================================
this_file = mfilename('fullpath');

if isempty(this_file)
    script_dir = pwd;
else
    script_dir = fileparts(this_file);
end

%% ========================================================================
% 1. FIGURE STYLE -- EDIT THESE VALUES IN THE FUTURE
% =========================================================================
%
% -------------------------------------------------------------------------
% FONT FAMILY
% -------------------------------------------------------------------------
STYLE.FONT_NAME = 'Times New Roman';   % <-- CHANGE ALL FIGURE FONTS HERE
%
% -------------------------------------------------------------------------
% FONT SIZES
% -------------------------------------------------------------------------
STYLE.TICK_FONT_SIZE   = 17;           % <-- tick-label font size
STYLE.AXIS_LABEL_SIZE  = 21;           % <-- x/y-axis label font size
STYLE.LEGEND_FONT_SIZE = 15;           % <-- legend font size
%
% -------------------------------------------------------------------------
% LINE / MARKER THICKNESS
% -------------------------------------------------------------------------
STYLE.LINE_WIDTH       = 3.2;          % <-- MAIN curve thickness
STYLE.ERROR_LINE_WIDTH = 2.3;          % <-- error-bar line thickness
STYLE.MARKER_SIZE      = 9.5;          % <-- marker size
STYLE.AXIS_LINE_WIDTH  = 1.8;          % <-- plot-box / axis thickness
STYLE.CAP_SIZE         = 10;           % <-- error-bar cap size
%
% -------------------------------------------------------------------------
% FIGURE SIZE / EXPORT QUALITY
% -------------------------------------------------------------------------
STYLE.FIGURE_POSITION = [100 100 1000 650];
STYLE.OUTPUT_DPI      = 400;           % <-- PNG export resolution
%
% -------------------------------------------------------------------------
% GRID APPEARANCE
% -------------------------------------------------------------------------
STYLE.GRID_ALPHA       = 0.16;
STYLE.MINOR_GRID_ALPHA = 0.08;
%
% -------------------------------------------------------------------------
% ETA TICK LABELS
% -------------------------------------------------------------------------
% We intentionally show only a sparse subset so the labels do NOT overlap.
STYLE.ETA_MAJOR_TICKS = [ ...
    0.3162277660, ...
    0.5623413252, ...
    1.0000000000, ...
    1.7782794100, ...
    2.2387211386];

%% ========================================================================
% 2. LOCATE THE THREE COMPLETED RESULT FOLDERS
% =========================================================================
DIR1 = locate_latest_result_folder( ...
    script_dir, ...
    'results_01_noise_tau_ablation_');

DIR2 = locate_latest_result_folder( ...
    script_dir, ...
    'results_02_q_tau_ablation_');

DIR3 = locate_latest_result_folder( ...
    script_dir, ...
    'results_03_sigma_tau_ablation_');

fprintf('\nUsing completed result folders:\n');
fprintf('  Noise sweep : %s\n',DIR1);
fprintf('  q sweep     : %s\n',DIR2);
fprintf('  sigma sweep : %s\n\n',DIR3);

%% ========================================================================
% 3. OUTPUT FOLDER
% =========================================================================
OUTPUT_DIR = fullfile(script_dir,'submission_figures_MoG-TeN');

if ~isfolder(OUTPUT_DIR)
    mkdir(OUTPUT_DIR);
end

%% ========================================================================
% 4. FIGURE 1 -- NOISE SWEEP
% =========================================================================
%
% Expected CSV from RUNME 1:
%   noise_sweep_summary.csv
%
T1_file = fullfile(DIR1,'noise_sweep_summary.csv');

assert(isfile(T1_file), ...
    'Could not find: %s',T1_file);

T1 = readtable(T1_file);
T1.Method = string(T1.Method);

fig = new_submission_figure(STYLE);
ax = axes(fig);
hold(ax,'on');

plot_eta_method( ...
    ax,T1,"proposed", ...
    'MoG-TeN', ...
    '-o',STYLE);

plot_eta_method( ...
    ax,T1,"tau0", ...
    'MoG-TeN ($\tau=0$)', ...
    '--s',STYLE);

plot_eta_method( ...
    ax,T1,"btd", ...
    'BTD', ...
    '-.^',STYLE);

plot_eta_method( ...
    ax,T1,"cp", ...
    'CPD', ...
    ':d',STYLE);

format_submission_axes(ax,STYLE);
format_eta_axis(ax,STYLE,T1.Eta);

set_axis_labels( ...
    ax, ...
    'Noise level, $\eta$', ...
    'Tensor $\mathcal{X}$ Reconstruction NMSE', ...
    STYLE);

% No title: the LaTeX subcaption will provide the experiment description.
make_submission_legend(ax,'northwest',STYLE);

save_submission_figure( ...
    fig, ...
    fullfile(OUTPUT_DIR,'submission_noise_sweep_X_NMSE'), ...
    STYLE);

close(fig);

%% ========================================================================
% 5. FIGURE 2 -- TRUE-RANK q SWEEP
% =========================================================================
%
% Expected CSV from RUNME 2:
%   all_q_summary_mean_std.csv
%
T2_file = fullfile(DIR2,'all_q_summary_mean_std.csv');

assert(isfile(T2_file), ...
    'Could not find: %s',T2_file);

T2 = readtable(T2_file);
T2.Method = string(T2.Method);

q_values = unique(T2.TrueRank,'sorted').';

fig = new_submission_figure(STYLE);
ax = axes(fig);
hold(ax,'on');

plot_q_method( ...
    ax,T2,"proposed", ...
    'MoG-TeN', ...
    '-o',STYLE);

plot_q_method( ...
    ax,T2,"tau0", ...
    'MoG-TeN ($\tau=0$)', ...
    '--s',STYLE);

plot_q_method( ...
    ax,T2,"btd", ...
    'BTD', ...
    '-.^',STYLE);

plot_q_method( ...
    ax,T2,"cp", ...
    'CPD', ...
    ':d',STYLE);

format_submission_axes(ax,STYLE);

xticks(ax,q_values);
xticklabels(ax,string(q_values));

if numel(q_values)>1
    qpad = 0.04*(max(q_values)-min(q_values));
    xlim(ax,[min(q_values)-qpad,max(q_values)+qpad]);
end

set_axis_labels( ...
    ax, ...
    'True template rank, $q$', ...
    'Tensor $\mathcal{X}$ Reconstruction NMSE', ...
    STYLE);

make_submission_legend(ax,'best',STYLE);

save_submission_figure( ...
    fig, ...
    fullfile(OUTPUT_DIR,'submission_q_sweep_X_NMSE'), ...
    STYLE);

close(fig);

%% ========================================================================
% 6. FIGURE 3 -- TEMPORAL-OVERLAP sigma SWEEP
% =========================================================================
%
% Expected CSV from RUNME 3:
%   summary_mean_std.csv
%
T3_file = fullfile(DIR3,'summary_mean_std.csv');

assert(isfile(T3_file), ...
    'Could not find: %s',T3_file);

T3 = readtable(T3_file);
T3.Method = string(T3.Method);

sigma_values = unique(T3.Sigma,'sorted').';

fig = new_submission_figure(STYLE);
ax = axes(fig);
hold(ax,'on');

plot_sigma_method( ...
    ax,T3,"proposed", ...
    'MoG-TeN', ...
    '-o',STYLE);

plot_sigma_method( ...
    ax,T3,"tau0", ...
    'MoG-TeN ($\tau=0$)', ...
    '--s',STYLE);

plot_sigma_method( ...
    ax,T3,"btd", ...
    'BTD', ...
    '-.^',STYLE);

plot_sigma_method( ...
    ax,T3,"cp", ...
    'CPD', ...
    ':d',STYLE);

format_submission_axes(ax,STYLE);

xticks(ax,sigma_values);
xticklabels(ax, ...
    arrayfun(@(x)sprintf('%.1f',x), ...
    sigma_values,'UniformOutput',false));

if numel(sigma_values)>1
    spad = 0.025*(max(sigma_values)-min(sigma_values));
    xlim(ax,[min(sigma_values)-spad,max(sigma_values)+spad]);
end

set_axis_labels( ...
    ax, ...
    'Overlap width, $\sigma$', ...
    'Tensor $\mathcal{X}$ Reconstruction NMSE', ...
    STYLE);

make_submission_legend(ax,'best',STYLE);

save_submission_figure( ...
    fig, ...
    fullfile(OUTPUT_DIR,'submission_sigma_sweep_X_NMSE'), ...
    STYLE);

close(fig);

%% ========================================================================
% 7. DONE
% =========================================================================
fprintf('\n============================================================\n');
fprintf('NEW MoG-TeN SUBMISSION FIGURES CREATED\n');
fprintf('============================================================\n');
fprintf('Output folder:\n%s\n\n',OUTPUT_DIR);

fprintf('Created:\n');
fprintf('  submission_noise_sweep_X_NMSE.png\n');
fprintf('  submission_q_sweep_X_NMSE.png\n');
fprintf('  submission_sigma_sweep_X_NMSE.png\n');

fprintf('\nEach figure is also saved as .fig and vector .pdf.\n');


%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function result_dir = locate_latest_result_folder(script_dir,prefix)
%LOCATE_LATEST_RESULT_FOLDER
% Find the newest result folder matching the specified experiment prefix.

D = dir(fullfile(script_dir,[prefix '*']));
D = D([D.isdir]);

if isempty(D)
    error( ...
        ['Could not locate a completed result folder matching:\n' ...
         '  %s*\n' ...
         'under:\n' ...
         '  %s'], ...
        prefix,script_dir);
end

% Folder timestamps are YYYYMMDD_HHMMSS, 
[~,idx] = sort({D.name});
D = D(idx);

result_dir = fullfile(D(end).folder,D(end).name);
end


function fig = new_submission_figure(STYLE)

fig = figure( ...
    'Visible','off', ...
    'Color','w', ...
    'Position',STYLE.FIGURE_POSITION);
end


function format_submission_axes(ax,STYLE)

grid(ax,'on');
box(ax,'on');

set(ax, ...
    'FontName',STYLE.FONT_NAME, ...
    'FontSize',STYLE.TICK_FONT_SIZE, ...
    'FontWeight','bold', ...
    'LineWidth',STYLE.AXIS_LINE_WIDTH, ...
    'TickDir','out', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',STYLE.GRID_ALPHA, ...
    'MinorGridAlpha',STYLE.MINOR_GRID_ALPHA);
end


function set_axis_labels(ax,xtext,ytext,STYLE)

xlabel(ax,xtext, ...
    'Interpreter','latex', ...
    'FontName',STYLE.FONT_NAME, ...
    'FontSize',STYLE.AXIS_LABEL_SIZE, ...
    'FontWeight','bold');

ylabel(ax,ytext, ...
    'Interpreter','latex', ...
    'FontName',STYLE.FONT_NAME, ...
    'FontSize',STYLE.AXIS_LABEL_SIZE, ...
    'FontWeight','bold');
end


function format_eta_axis(ax,STYLE,eta_values)

eta_values = sort(unique(eta_values(:)));

set(ax,'XScale','log');

requested = STYLE.ETA_MAJOR_TICKS(:);

ticks = [];

for k = 1:numel(requested)

    [distance,idx] = ...
        min(abs(log(eta_values)-log(requested(k))));

    if distance<0.08
        ticks(end+1,1) = eta_values(idx); %#ok<AGROW>
    end
end

ticks = unique(ticks,'stable');

if numel(ticks)<3

    idx = unique(round(linspace( ...
        1,numel(eta_values), ...
        min(5,numel(eta_values)))));

    ticks = eta_values(idx);
end

xticks(ax,ticks);

xticklabels(ax, ...
    arrayfun(@(x)sprintf('%.3g',x), ...
    ticks,'UniformOutput',false));

xlim(ax,[ ...
    min(eta_values)/1.07, ...
    max(eta_values)*1.07]);
end


function lgd = make_submission_legend(ax,location,STYLE)

lgd = legend(ax, ...
    'Location',location, ...
    'Interpreter','latex', ...
    'FontName',STYLE.FONT_NAME, ...
    'FontSize',STYLE.LEGEND_FONT_SIZE, ...
    'Box','off');


try
    lgd.ItemTokenSize = [24 18];
catch
end
end


function plot_eta_method(ax,T,method,label,line_style,STYLE)

A = T(T.Method==string(method),:);
A = sortrows(A,'Eta');

plot(ax, ...
    A.Eta, ...
    A.CleanReconNMSE, ...
    line_style, ...
    'LineWidth',STYLE.LINE_WIDTH, ...
    'MarkerSize',STYLE.MARKER_SIZE, ...
    'DisplayName',label);
end


function plot_q_method(ax,T,method,label,line_style,STYLE)

A = T(T.Method==string(method),:);
A = sortrows(A,'TrueRank');

if ismember('CleanReconNMSE_Std',A.Properties.VariableNames)

    errorbar(ax, ...
        A.TrueRank, ...
        A.CleanReconNMSE_Mean, ...
        A.CleanReconNMSE_Std, ...
        line_style, ...
        'LineWidth',STYLE.LINE_WIDTH, ...
        'MarkerSize',STYLE.MARKER_SIZE, ...
        'CapSize',STYLE.CAP_SIZE, ...
        'DisplayName',label);

else

    plot(ax, ...
        A.TrueRank, ...
        A.CleanReconNMSE_Mean, ...
        line_style, ...
        'LineWidth',STYLE.LINE_WIDTH, ...
        'MarkerSize',STYLE.MARKER_SIZE, ...
        'DisplayName',label);
end
end


function plot_sigma_method(ax,T,method,label,line_style,STYLE)

A = T(T.Method==string(method),:);
A = sortrows(A,'Sigma');

if ismember('CleanReconNMSE_Std',A.Properties.VariableNames)

    errorbar(ax, ...
        A.Sigma, ...
        A.CleanReconNMSE_Mean, ...
        A.CleanReconNMSE_Std, ...
        line_style, ...
        'LineWidth',STYLE.LINE_WIDTH, ...
        'MarkerSize',STYLE.MARKER_SIZE, ...
        'CapSize',STYLE.CAP_SIZE, ...
        'DisplayName',label);

else

    plot(ax, ...
        A.Sigma, ...
        A.CleanReconNMSE_Mean, ...
        line_style, ...
        'LineWidth',STYLE.LINE_WIDTH, ...
        'MarkerSize',STYLE.MARKER_SIZE, ...
        'DisplayName',label);
end
end


function save_submission_figure(fig,basepath,STYLE)

% High-resolution PNG 
try
    exportgraphics( ...
        fig, ...
        [basepath '.png'], ...
        'Resolution',STYLE.OUTPUT_DPI);
catch
    print( ...
        fig, ...
        [basepath '.png'], ...
        '-dpng', ...
        sprintf('-r%d',STYLE.OUTPUT_DPI));
end

% Editable MATLAB figure.
try
    savefig(fig,[basepath '.fig']);
catch
end

% Vector PDF for publication if desired.
try
    exportgraphics( ...
        fig, ...
        [basepath '.pdf'], ...
        'ContentType','vector');
catch
    try
        print(fig,[basepath '.pdf'],'-dpdf','-painters');
    catch
    end
end
end
