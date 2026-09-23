%% RUNME_3_SIGMA_SWEEP_WITH_TAU_ABLATION.M
%
% EXPERIMENT 3 -- TEMPORAL-OVERLAP SIGMA SWEEP
% ================================================================
%
% RESWEEP ONLY.
%
% sigma = [1.5 2 2.5 3 3.5 4 4.5 5]
% eta   = 0.5623413252
% q     = 2
%
% METHODS
%   Proposed
%   Proposed (tau=0), paired nuclear-norm ablation
%   BTD
%   CPD
%
% At each sigma the positive-tau OURS model is independently reswept and
% selected by the balanced CLEAN component-correlation criterion.  Then only
% tau is changed to zero for the paired ablation.
%
% PRIMARY SUBMISSION FIGURE
%   submission_sigma_sweep_X_NMSE.png
%
% Additional C_r / alpha_r correlation diagnostics are saved separately.

% DONT FORGET TO REPLACE THE TENSORLAB PATH WITH THE USER'S OWN PATH 


clear; clc; close all;

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file), script_dir = pwd; else, script_dir = fileparts(this_file); end
addpath(genpath(fullfile(script_dir,'lib')));

CFG.tensorlabPath = ' '; % here put your tensorlab path 
assert(isfolder(CFG.tensorlabPath),'Tensorlab folder not found: %s',CFG.tensorlabPath);
addpath(genpath(CFG.tensorlabPath));
assert(exist('ll1','file')==2,'Tensorlab ll1.m not visible.');
assert(exist('cpd','file')==2,'Tensorlab cpd.m not visible.');

%% Data
DATA.N = 128;
DATA.T = 20;
DATA.R = 3;
DATA.q = 2;
DATA.target_C_fro = 10.0;
DATA.background_scale = 0.08;
DATA.seed_structure = 20262014;

TEMP.centers = [5.0 10.5 16.0];
OVERLAP.sigma_values = [1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0];

NOISE.eta = 0.5623413252;
NOISE.seed_resweep_base = 20264000;
NOISE.seed_eval_base = 20265000;

EVAL.n_repetitions = 3;

%% OURS search
SEARCH.n_random_design = 6;
SEARCH.seed_design = 20260915;
SEARCH.log10_gamma  = [-4, log10(0.5)];
SEARCH.log10_tau    = [-3, log10(2.0)];
SEARCH.log10_mu     = [-2, log10(50)];
SEARCH.log10_lambda = [-3, 0];

SOLVER.zeta = 2.0;
SOLVER.max_iter = 220;
SOLVER.min_iter = 15;
SOLVER.tol_primal = 1e-4;
SOLVER.tol_dual = 2e-4;
SOLVER.tol_obj = 1e-6;
SOLVER.alpha_tol = 1e-10;
SOLVER.range_tol = 1e-10;
SOLVER.eig_tol = 1e-12;
SOLVER.bisect_tol = 1e-12;
SOLVER.bisect_maxit = 100;
SOLVER.init_density = 0.15;
SOLVER.verbose_every = 0;

TUNE.n_restarts = 3;
TUNE.seed_solver0 = 97501;

OURS_EVAL.n_restarts = 3;
OURS_EVAL.seed0 = 98001;

%% BTD / CPD
BTD.n_restarts = 2;
BTD.seed0 = 99001;
BTD.max_iter = 350;
BTD.tol_fun = 1e-10;
BTD.tol_x = 1e-8;
BTD.display = 0;

CP.n_restarts = 2;
CP.seed0 = 99501;
CP.max_iter = 350;
CP.tol_fun = 1e-10;
CP.tol_x = 1e-8;
CP.display = 0;

%% Output
timestamp = datestr(now,'yyyymmdd_HHMMSS');
root_dir = fullfile(script_dir,['results_03_sigma_tau_ablation_' timestamp]);
fig_dir = fullfile(root_dir,'submission');
diag_dir = fullfile(root_dir,'diagnostics');
tuning_dir = fullfile(root_dir,'tuning');
submission_dir = fullfile(script_dir,'submission_figures');
mkdir(root_dir); mkdir(fig_dir); mkdir(diag_dir); mkdir(tuning_dir);
if ~isfolder(submission_dir), mkdir(submission_dir); end

try
    copyfile(this_file,fullfile(root_dir,'RUNME_3_sigma_sweep_with_tau_ablation.m'));
catch
end

%% Shared clean structural templates
blocks = make_component_blocks(DATA.N,DATA.R);
C_raw = generate_latent_group_templates( ...
    DATA.N,DATA.R,DATA.q,blocks, ...
    DATA.background_scale,DATA.target_C_fro,DATA.seed_structure);

configs = build_joint_hyperparameter_design(SEARCH,SOLVER.zeta);
writetable(configs,fullfile(root_dir,'positive_tau_hyperparameter_candidates.csv'));

[L,V_L,ell_L] = temporal_laplacian(DATA.T);

selected_full_hyp = cell(numel(OVERLAP.sigma_values),1);
selected_tau0_hyp = cell(numel(OVERLAP.sigma_values),1);
selected_hyp_rows = table();

%% Resweep tuning at every sigma
for sidx = 1:numel(OVERLAP.sigma_values)

    sigma = OVERLAP.sigma_values(sidx);

    alpha_raw = make_gaussian_temporal_factors( ...
        DATA.T,DATA.R,TEMP.centers,sigma);

    Xclean = reconstruct_tensor(C_raw,alpha_raw);
    [Ctrue,atrue] = canonicalize_truth(C_raw,alpha_raw);

    Htune = generate_template_noise_directions( ...
        DATA.N,DATA.R,NOISE.seed_resweep_base+sidx);

    [Xtune,~,~,~,~] = add_template_noise_at_eta( ...
        C_raw,alpha_raw,Htune,NOISE.eta);

    fprintf('\n============================================================\n');
    fprintf('SIGMA RESWEEP | sigma=%.2f\n',sigma);
    fprintf('============================================================\n');

    [~,~,best_corr,hyp_full,cfg_table,restart_table] = ...
        sweep_ours_hyperparameters_by_correlation( ...
        Xtune,Xclean,Ctrue,atrue, ...
        DATA.R,L,V_L,ell_L,configs,SOLVER,TUNE.n_restarts, ...
        TUNE.seed_solver0+100*sidx);

    hyp_tau0 = hyp_full;
    hyp_tau0.tau = 0;

    selected_full_hyp{sidx} = hyp_full;
    selected_tau0_hyp{sidx} = hyp_tau0;

    sdir = fullfile(tuning_dir,sprintf('sigma_%0.2f',sigma));
    mkdir(sdir);
    writetable(cfg_table,fullfile(sdir,'positive_tau_config_search.csv'));
    writetable(restart_table,fullfile(sdir,'positive_tau_all_restarts.csv'));

    selected_hyp_rows = [selected_hyp_rows; table( ...
        sigma,hyp_full.ConfigID,hyp_full.gamma,hyp_full.tau,hyp_full.mu, ...
        hyp_full.lambda,hyp_full.zeta,best_corr.C_mean,best_corr.alpha_mean, ...
        best_corr.balanced_score, ...
        'VariableNames',{'Sigma','ConfigID','gamma','tau','mu','lambda','zeta', ...
        'MeanCCorrelation','MeanAlphaCorrelation','BalancedCorrelationScore'})]; %#ok<AGROW>

    writetable(selected_hyp_rows,fullfile(root_dir,'selected_hyperparameters_partial.csv'));
end

writetable(selected_hyp_rows,fullfile(root_dir,'selected_hyperparameters_by_sigma.csv'));

%% Held-out comparison
raw_rows = table();
component_rows = table();
Lvec = DATA.q*ones(1,DATA.R);
Pcp = DATA.R*DATA.q;

for sidx = 1:numel(OVERLAP.sigma_values)

    sigma = OVERLAP.sigma_values(sidx);

    alpha_raw = make_gaussian_temporal_factors( ...
        DATA.T,DATA.R,TEMP.centers,sigma);

    Xclean = reconstruct_tensor(C_raw,alpha_raw);
    [Ctrue,atrue] = canonicalize_truth(C_raw,alpha_raw);

    hyp_full = selected_full_hyp{sidx};
    hyp_tau0 = selected_tau0_hyp{sidx};

    for rep = 1:EVAL.n_repetitions

        Heval = generate_template_noise_directions( ...
            DATA.N,DATA.R,NOISE.seed_eval_base+rep);

        [Xobs,~,~,actual_eta,~] = add_template_noise_at_eta( ...
            C_raw,alpha_raw,Heval,NOISE.eta);

        seed_eval = OURS_EVAL.seed0 + 10000*sidx + 100*rep;

        [ours_result,~] = fit_ours_multistart( ...
            Xobs,DATA.R,L,V_L,ell_L,hyp_full,SOLVER, ...
            OURS_EVAL.n_restarts,seed_eval);

        [tau0_result,~] = fit_ours_multistart( ...
            Xobs,DATA.R,L,V_L,ell_L,hyp_tau0,SOLVER, ...
            OURS_EVAL.n_restarts,seed_eval);

        ours_metrics = evaluate_solution_general( ...
            Xobs,Xclean,ours_result.C,ours_result.alpha,Ctrue,atrue,false);
        ours_corr = factor_correlations_balanced( ...
            Ctrue,atrue,ours_result.C,ours_result.alpha,false);

        tau0_metrics = evaluate_solution_general( ...
            Xobs,Xclean,tau0_result.C,tau0_result.alpha,Ctrue,atrue,false);
        tau0_corr = factor_correlations_balanced( ...
            Ctrue,atrue,tau0_result.C,tau0_result.alpha,false);

        btd_result = fit_btd_multistart( ...
            Xobs,Lvec,BTD,BTD.seed0+10000*sidx+100*rep);
        btd_metrics = evaluate_solution_general( ...
            Xobs,Xclean,btd_result.C,btd_result.alpha,Ctrue,atrue,true);
        btd_corr = factor_correlations_balanced( ...
            Ctrue,atrue,btd_result.C,btd_result.alpha,true);

        cp_result = fit_cp_multistart( ...
            Xobs,Pcp,CP,CP.seed0+10000*sidx+100*rep);
        cp_metrics = evaluate_tensor_reconstruction_only(Xobs,Xclean,cp_result.Xhat);

        raw_rows = append_summary(raw_rows,sigma,rep,actual_eta,'proposed', ...
            ours_metrics.clean_recon_nmse,ours_corr.C_mean, ...
            ours_corr.alpha_mean,ours_corr.balanced_score);

        raw_rows = append_summary(raw_rows,sigma,rep,actual_eta,'tau0', ...
            tau0_metrics.clean_recon_nmse,tau0_corr.C_mean, ...
            tau0_corr.alpha_mean,tau0_corr.balanced_score);

        raw_rows = append_summary(raw_rows,sigma,rep,actual_eta,'btd', ...
            btd_metrics.clean_recon_nmse,btd_corr.C_mean, ...
            btd_corr.alpha_mean,btd_corr.balanced_score);

        raw_rows = append_summary(raw_rows,sigma,rep,actual_eta,'cp', ...
            cp_metrics.clean_recon_nmse,NaN,NaN,NaN);

        for r = 1:DATA.R
            component_rows = append_component(component_rows,sigma,rep,'proposed',r, ...
                ours_corr.C_each(r),ours_corr.alpha_each(r));
            component_rows = append_component(component_rows,sigma,rep,'tau0',r, ...
                tau0_corr.C_each(r),tau0_corr.alpha_each(r));
            component_rows = append_component(component_rows,sigma,rep,'btd',r, ...
                btd_corr.C_each(r),btd_corr.alpha_each(r));
        end
    end

    writetable(raw_rows,fullfile(root_dir,'heldout_summary_partial.csv'));
    writetable(component_rows,fullfile(root_dir,'component_correlations_partial.csv'));
end

writetable(raw_rows,fullfile(root_dir,'heldout_summary.csv'));
writetable(component_rows,fullfile(root_dir,'component_correlations.csv'));

S = aggregate_by_sigma(raw_rows,OVERLAP.sigma_values);
writetable(S,fullfile(root_dir,'summary_mean_std.csv'));

%% Submission figure
fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

plot_sigma_method(ax,S,"proposed",'Proposed','-o');
plot_sigma_method(ax,S,"tau0",'Proposed ($\tau=0$)','--s');
plot_sigma_method(ax,S,"btd",'BTD','-.^');
plot_sigma_method(ax,S,"cp",'CPD',':d');

submission_axes(ax);
xticks(ax,OVERLAP.sigma_values);
xticklabels(ax,arrayfun(@(x)sprintf('%.1f',x),OVERLAP.sigma_values,'UniformOutput',false));
xlim(ax,[1.4 5.1]);

submission_labels(ax,'Overlap width, $\sigma$','Tensor $\mathcal{X}$ Reconstruction NMSE');
submission_legend(ax,'best');

% Intentionally no title.
save_figure_pair(fig,fullfile(fig_dir,'submission_sigma_sweep_X_NMSE'));
save_figure_pair(fig,fullfile(submission_dir,'submission_sigma_sweep_X_NMSE'));
close(fig);

%% Diagnostics
plot_sigma_corr(S,OVERLAP.sigma_values,'C', ...
    fullfile(diag_dir,'diagnostic_C_correlation_vs_sigma'));

plot_sigma_corr(S,OVERLAP.sigma_values,'alpha', ...
    fullfile(diag_dir,'diagnostic_alpha_correlation_vs_sigma'));

plot_sigma_balanced(S,OVERLAP.sigma_values, ...
    fullfile(diag_dir,'diagnostic_balanced_factor_score_vs_sigma'));

fprintf('\nTASK 3 COMPLETE\nResults: %s\n',root_dir);


%% ========================================================================
% LOCAL HELPERS
% =========================================================================

function T = append_summary(T,sigma,rep,eta,method,clean_nmse,ccorr,acorr,score)
row = table(sigma,rep,eta,string(method),clean_nmse,ccorr,acorr,score, ...
    'VariableNames',{'Sigma','Repetition','Eta','Method','CleanReconNMSE', ...
    'MeanCCorrelation','MeanAlphaCorrelation','BalancedCorrelationScore'});
T = [T;row];
end


function T = append_component(T,sigma,rep,method,r,ccorr,acorr)
row = table(sigma,rep,string(method),r,ccorr,acorr, ...
    'VariableNames',{'Sigma','Repetition','Method','Component', ...
    'CCorrelation','AlphaCorrelation'});
T = [T;row];
end


function S = aggregate_by_sigma(T,sigma_values)
methods = unique(T.Method,'stable');
S = table();

for m = 1:numel(methods)
    for sigma = sigma_values
        A = T(T.Method==methods(m) & abs(T.Sigma-sigma)<1e-12,:);
        row = table(methods(m),sigma, ...
            mean(A.CleanReconNMSE,'omitnan'),std(A.CleanReconNMSE,'omitnan'), ...
            mean(A.MeanCCorrelation,'omitnan'),std(A.MeanCCorrelation,'omitnan'), ...
            mean(A.MeanAlphaCorrelation,'omitnan'),std(A.MeanAlphaCorrelation,'omitnan'), ...
            mean(A.BalancedCorrelationScore,'omitnan'),std(A.BalancedCorrelationScore,'omitnan'), ...
            'VariableNames',{'Method','Sigma','CleanReconNMSE_Mean', ...
            'CleanReconNMSE_Std','CCorrelation_Mean','CCorrelation_Std', ...
            'AlphaCorrelation_Mean','AlphaCorrelation_Std', ...
            'BalancedScore_Mean','BalancedScore_Std'});
        S = [S;row];
    end
end
end


function plot_sigma_method(ax,S,method,label,style)
A = sortrows(S(S.Method==string(method),:),'Sigma');
errorbar(ax,A.Sigma,A.CleanReconNMSE_Mean,A.CleanReconNMSE_Std,style, ...
    'LineWidth',2.2,'MarkerSize',7,'CapSize',7,'DisplayName',label);
end


function plot_sigma_corr(S,sigma_values,which_factor,basepath)

methods = ["proposed","tau0","btd"];
names = {'Proposed','Proposed ($\tau=0$)','BTD'};
styles = {'-o','--s','-.^'};

fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

for m = 1:numel(methods)
    A = sortrows(S(S.Method==methods(m),:),'Sigma');

    if strcmp(which_factor,'C')
        yy = A.CCorrelation_Mean; ee = A.CCorrelation_Std;
    else
        yy = A.AlphaCorrelation_Mean; ee = A.AlphaCorrelation_Std;
    end

    errorbar(ax,A.Sigma,yy,ee,styles{m}, ...
        'LineWidth',1.8,'MarkerSize',7,'CapSize',7,'DisplayName',names{m});
end

submission_axes(ax);
xticks(ax,sigma_values);
xticklabels(ax,arrayfun(@(x)sprintf('%.1f',x),sigma_values,'UniformOutput',false));

if strcmp(which_factor,'C')
    submission_labels(ax,'Overlap width, $\sigma$','Mean $C_r$ Correlation');
else
    submission_labels(ax,'Overlap width, $\sigma$','Mean $\alpha_r$ Correlation');
end

ylim(ax,[-0.05 1.05]);
submission_legend(ax,'best');
save_figure_pair(fig,basepath);
close(fig);
end


function plot_sigma_balanced(S,sigma_values,basepath)

methods = ["proposed","tau0","btd"];
names = {'Proposed','Proposed ($\tau=0$)','BTD'};
styles = {'-o','--s','-.^'};

fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

for m = 1:numel(methods)
    A = sortrows(S(S.Method==methods(m),:),'Sigma');
    errorbar(ax,A.Sigma,A.BalancedScore_Mean,A.BalancedScore_Std,styles{m}, ...
        'LineWidth',1.8,'MarkerSize',7,'CapSize',7,'DisplayName',names{m});
end

submission_axes(ax);
xticks(ax,sigma_values);
xticklabels(ax,arrayfun(@(x)sprintf('%.1f',x),sigma_values,'UniformOutput',false));
submission_labels(ax,'Overlap width, $\sigma$','$\min(\bar{\rho}_C,\bar{\rho}_{\alpha})$');
ylim(ax,[-0.05 1.05]);
submission_legend(ax,'best');
save_figure_pair(fig,basepath);
close(fig);
end
