%% RUNME_1_NOISE_SWEEP_WITH_TAU_ABLATION.M
%
% EXPERIMENT 1 -- TEMPLATE-NOISE SWEEP, q=2
% ================================================================
%
% METHODS
%   Proposed                         : selected positive tau
%   Proposed (tau=0)                 : paired nuclear-norm ablation
%   BTD                              : L=[2 2 2]
%   CPD                              : P=6
%
% OURS HYPERPARAMETER SELECTION
% -----------------------------
% At EACH eta, jointly search gamma/tau/mu/lambda for the positive-tau
% model.  The selected tuple maximizes
%
%   min( mean_r corr(C_r^clean,C_r^estimated),
%        mean_r corr(alpha_r^clean,alpha_r^estimated) ).
%
% Clean-X reconstruction is NOT used for hyperparameter selection.
%
% After the tuple is selected, create the paired ablation by setting ONLY
%
%   tau -> 0
%
% while keeping gamma, mu, lambda and zeta fixed.
%
% PRIMARY SUBMISSION FIGURE
% -------------------------
% submission_noise_sweep_X_NMSE.png
%
% It has NO title and compares:
%   Proposed / Proposed (tau=0) / BTD / CPD
%
% Additional C_r / alpha_r correlation figures are saved under diagnostics/.
%
% DATA GENERATION IS THE SAME AS THE CURRENT PAPER SUITE.
%
% RUN:
%   RUNME_1_noise_sweep_with_tau_ablation

% DONT FORGET TO REPLACE THE TENSORLAB PATH WITH THE USER'S OWN PATH 


clear; clc; close all;

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file), script_dir = pwd; else, script_dir = fileparts(this_file); end
addpath(genpath(fullfile(script_dir,'lib')));

CFG.tensorlabPath = ' ';  % to be replaced with the user's own tensorlab path 
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
DATA.seed_structure = 20262002;

TEMP.centers = [5.0 10.5 16.0];
TEMP.sigma = 3.0;

ETA.values = [ ...
    0.3162277660, ...
    0.4216965034, ...
    0.5623413252, ...
    0.7498942093, ...
    1.0000000000, ...
    1.3335214322, ...
    1.7782794100, ...
    2.2387211386];

ETA.equivalent_snr_db = [10 7.5 5 2.5 0 -2.5 -5 -7];
NOISE.seed_evaluation = 20263002;

%% OURS search
SEARCH.n_random_design = 6;
SEARCH.seed_design = 20260915;
SEARCH.log10_gamma  = [-4, log10(0.5)];
SEARCH.log10_tau    = [-3, log10(2.0)];
SEARCH.log10_mu     = [-2, log10(50)];
SEARCH.log10_lambda = [-3, 0];

EVAL.n_restarts = 3;
EVAL.seed_solver0 = 98501;

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
root_dir = fullfile(script_dir,['results_01_noise_tau_ablation_' timestamp]);
fig_dir = fullfile(root_dir,'submission');
diag_dir = fullfile(root_dir,'diagnostics');
submission_dir = fullfile(script_dir,'submission_figures');
mkdir(root_dir); mkdir(fig_dir); mkdir(diag_dir);
if ~isfolder(submission_dir), mkdir(submission_dir); end

try
    copyfile(this_file,fullfile(root_dir,'RUNME_1_noise_sweep_with_tau_ablation.m'));
catch
end

%% Clean truth
blocks = make_component_blocks(DATA.N,DATA.R);
C_true_raw = generate_latent_group_templates( ...
    DATA.N,DATA.R,DATA.q,blocks, ...
    DATA.background_scale,DATA.target_C_fro,DATA.seed_structure);

alpha_true_raw = make_gaussian_temporal_factors( ...
    DATA.T,DATA.R,TEMP.centers,TEMP.sigma);

X_clean = reconstruct_tensor(C_true_raw,alpha_true_raw);
[C_true,alpha_true,alpha_l2_before,canon_error] = ...
    canonicalize_truth(C_true_raw,alpha_true_raw);

assert(canon_error<1e-10,'Truth canonicalization error is unexpectedly large.');

save(fullfile(root_dir,'clean_ground_truth.mat'), ...
    'C_true_raw','C_true','alpha_true_raw','alpha_true', ...
    'alpha_l2_before','canon_error','X_clean','-v7.3');

%% Shared noise direction and hyperparameter candidates
H_eval = generate_template_noise_directions(DATA.N,DATA.R,NOISE.seed_evaluation);
configs = build_joint_hyperparameter_design(SEARCH,SOLVER.zeta);
writetable(configs,fullfile(root_dir,'positive_tau_hyperparameter_candidates.csv'));

[L,V_L,ell_L] = temporal_laplacian(DATA.T);
Lvec = DATA.q*ones(1,DATA.R);
Pcp = DATA.R*DATA.q;

summary_rows = table();
component_rows = table();
selected_hyp_rows = table();

%% Eta sweep
for eidx = 1:numel(ETA.values)

    eta = ETA.values(eidx);
    snr_db = ETA.equivalent_snr_db(eidx);

    eta_dir = fullfile(root_dir,eta_folder_name(eta));
    mkdir(eta_dir);

    [Xobs,C_corrupt_raw,NoiseTensor,actual_eta,actual_snr_db] = ...
        add_template_noise_at_eta(C_true_raw,alpha_true_raw,H_eval,eta);

    save(fullfile(eta_dir,'generated_observation.mat'), ...
        'Xobs','X_clean','C_corrupt_raw','NoiseTensor', ...
        'eta','actual_eta','actual_snr_db','snr_db','-v7.3');

    fprintf('\n============================================================\n');
    fprintf('NOISE SWEEP | eta=%.6g\n',eta);
    fprintf('============================================================\n');

    seed_ours = EVAL.seed_solver0 + 100*eidx;

    %% Proposed: positive tau, selected by component correlations
    [ours_result,ours_metrics,ours_corr,ours_hyp, ...
        ours_config_table,ours_restart_table] = ...
        sweep_ours_hyperparameters_by_correlation( ...
        Xobs,X_clean,C_true,alpha_true, ...
        DATA.R,L,V_L,ell_L,configs,SOLVER,EVAL.n_restarts,seed_ours);

    writetable(ours_config_table,fullfile(eta_dir,'positive_tau_config_search.csv'));
    writetable(ours_restart_table,fullfile(eta_dir,'positive_tau_all_restarts.csv'));

    %% Paired nuclear ablation: ONLY tau -> 0
    tau0_hyp = ours_hyp;
    tau0_hyp.tau = 0;

    [tau0_result,tau0_runs] = fit_ours_multistart( ...
        Xobs,DATA.R,L,V_L,ell_L,tau0_hyp,SOLVER,EVAL.n_restarts,seed_ours);

    tau0_metrics = evaluate_solution_general( ...
        Xobs,X_clean,tau0_result.C,tau0_result.alpha,C_true,alpha_true,false);

    tau0_corr = factor_correlations_balanced( ...
        C_true,alpha_true,tau0_result.C,tau0_result.alpha,false);

    writetable(tau0_runs,fullfile(eta_dir,'tau0_paired_restarts.csv'));

    %% BTD
    btd_result = fit_btd_multistart( ...
        Xobs,Lvec,BTD,BTD.seed0+100*eidx);

    btd_metrics = evaluate_solution_general( ...
        Xobs,X_clean,btd_result.C,btd_result.alpha,C_true,alpha_true,true);

    btd_corr = factor_correlations_balanced( ...
        C_true,alpha_true,btd_result.C,btd_result.alpha,true);

    %% CPD
    cp_result = fit_cp_multistart( ...
        Xobs,Pcp,CP,CP.seed0+100*eidx);

    cp_metrics = evaluate_tensor_reconstruction_only(Xobs,X_clean,cp_result.Xhat);

    %% Save full selected solutions
    save(fullfile(eta_dir,'selected_solutions.mat'), ...
        'ours_result','ours_metrics','ours_corr','ours_hyp', ...
        'tau0_result','tau0_metrics','tau0_corr','tau0_hyp', ...
        'btd_result','btd_metrics','btd_corr', ...
        'cp_result','cp_metrics','C_true','alpha_true','X_clean','Xobs','-v7.3');

    %% Summary
    summary_rows = append_summary(summary_rows,eta,snr_db,'proposed', ...
        ours_metrics.clean_recon_nmse,ours_metrics.noisy_recon_nmse, ...
        ours_corr.C_mean,ours_corr.alpha_mean,ours_corr.balanced_score);

    summary_rows = append_summary(summary_rows,eta,snr_db,'tau0', ...
        tau0_metrics.clean_recon_nmse,tau0_metrics.noisy_recon_nmse, ...
        tau0_corr.C_mean,tau0_corr.alpha_mean,tau0_corr.balanced_score);

    summary_rows = append_summary(summary_rows,eta,snr_db,'btd', ...
        btd_metrics.clean_recon_nmse,btd_metrics.noisy_recon_nmse, ...
        btd_corr.C_mean,btd_corr.alpha_mean,btd_corr.balanced_score);

    summary_rows = append_summary(summary_rows,eta,snr_db,'cp', ...
        cp_metrics.clean_recon_nmse,cp_metrics.noisy_recon_nmse,NaN,NaN,NaN);

    for r = 1:DATA.R
        component_rows = append_component(component_rows,eta,'proposed',r, ...
            ours_corr.C_each(r),ours_corr.alpha_each(r));
        component_rows = append_component(component_rows,eta,'tau0',r, ...
            tau0_corr.C_each(r),tau0_corr.alpha_each(r));
        component_rows = append_component(component_rows,eta,'btd',r, ...
            btd_corr.C_each(r),btd_corr.alpha_each(r));
    end

    selected_hyp_rows = [selected_hyp_rows; table( ...
        eta,ours_hyp.ConfigID,ours_hyp.gamma,ours_hyp.tau,ours_hyp.mu, ...
        ours_hyp.lambda,ours_hyp.zeta,ours_corr.C_mean,ours_corr.alpha_mean, ...
        ours_corr.balanced_score, ...
        'VariableNames',{'Eta','ConfigID','gamma','tau','mu','lambda','zeta', ...
        'MeanCCorrelation','MeanAlphaCorrelation','BalancedCorrelationScore'})]; %#ok<AGROW>

    writetable(summary_rows,fullfile(root_dir,'summary_partial.csv'));
    writetable(component_rows,fullfile(root_dir,'component_correlations_partial.csv'));
    writetable(selected_hyp_rows,fullfile(root_dir,'selected_hyperparameters_partial.csv'));
end

%% Final CSVs
writetable(summary_rows,fullfile(root_dir,'noise_sweep_summary.csv'));
writetable(component_rows,fullfile(root_dir,'noise_sweep_component_correlations.csv'));
writetable(selected_hyp_rows,fullfile(root_dir,'selected_hyperparameters_by_eta.csv'));

%% Submission figure: X_clean reconstruction
fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

plot_method_eta(ax,summary_rows,"proposed",'Proposed','-o');
plot_method_eta(ax,summary_rows,"tau0",'Proposed ($\tau=0$)','--s');
plot_method_eta(ax,summary_rows,"btd",'BTD','-.^');
plot_method_eta(ax,summary_rows,"cp",'CPD',':d');

submission_axes(ax);
paper_eta_axis(ax,ETA.values);
submission_labels(ax,'Noise level, $\eta$','Tensor $\mathcal{X}$ Reconstruction NMSE');
submission_legend(ax,'northwest');

% Intentionally no title.
save_figure_pair(fig,fullfile(fig_dir,'submission_noise_sweep_X_NMSE'));
save_figure_pair(fig,fullfile(submission_dir,'submission_noise_sweep_X_NMSE'));
close(fig);

%% Diagnostics: component correlations
plot_corr_diagnostic(component_rows,ETA.values,'C', ...
    fullfile(diag_dir,'diagnostic_C_correlation_vs_eta'));

plot_corr_diagnostic(component_rows,ETA.values,'alpha', ...
    fullfile(diag_dir,'diagnostic_alpha_correlation_vs_eta'));

plot_balanced_diagnostic(summary_rows,ETA.values, ...
    fullfile(diag_dir,'diagnostic_balanced_factor_score_vs_eta'));

fprintf('\nTASK 1 COMPLETE\nResults: %s\n',root_dir);


%% ========================================================================
% LOCAL HELPERS
% =========================================================================

function T = append_summary(T,eta,snr_db,method,clean_nmse,observed_nmse,ccorr,acorr,score)
row = table(eta,snr_db,string(method),clean_nmse,observed_nmse,ccorr,acorr,score, ...
    'VariableNames',{'Eta','EquivalentSNR_dB','Method','CleanReconNMSE', ...
    'ObservedReconNMSE','MeanCCorrelation','MeanAlphaCorrelation', ...
    'BalancedCorrelationScore'});
T = [T;row];
end


function T = append_component(T,eta,method,r,ccorr,acorr)
row = table(eta,string(method),r,ccorr,acorr, ...
    'VariableNames',{'Eta','Method','Component','CCorrelation','AlphaCorrelation'});
T = [T;row];
end


function plot_method_eta(ax,T,method,label,style)
A = sortrows(T(T.Method==string(method),:),'Eta');
plot(ax,A.Eta,A.CleanReconNMSE,style, ...
    'LineWidth',2.2,'MarkerSize',7,'DisplayName',label);
end


function plot_corr_diagnostic(T,eta_values,which_factor,basepath)

methods = ["proposed","tau0","btd"];
names = {'Proposed','Proposed ($\tau=0$)','BTD'};
styles = {'-o','--s','-.^'};

fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

for m = 1:numel(methods)

    means = nan(size(eta_values));
    stds = nan(size(eta_values));

    for k = 1:numel(eta_values)
        A = T(T.Method==methods(m) & abs(T.Eta-eta_values(k))<1e-10,:);
        if strcmp(which_factor,'C')
            vals = A.CCorrelation;
        else
            vals = A.AlphaCorrelation;
        end
        means(k) = mean(vals,'omitnan');
        stds(k) = std(vals,'omitnan');
    end

    errorbar(ax,eta_values,means,stds,styles{m}, ...
        'LineWidth',1.8,'MarkerSize',7,'DisplayName',names{m});
end

submission_axes(ax);
paper_eta_axis(ax,eta_values);

if strcmp(which_factor,'C')
    submission_labels(ax,'Noise level, $\eta$','Mean $C_r$ Correlation');
else
    submission_labels(ax,'Noise level, $\eta$','Mean $\alpha_r$ Correlation');
end

ylim(ax,[-0.05 1.05]);
submission_legend(ax,'best');
save_figure_pair(fig,basepath);
close(fig);
end


function plot_balanced_diagnostic(T,eta_values,basepath)

methods = ["proposed","tau0","btd"];
names = {'Proposed','Proposed ($\tau=0$)','BTD'};
styles = {'-o','--s','-.^'};

fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

for m = 1:numel(methods)
    A = sortrows(T(T.Method==methods(m),:),'Eta');
    plot(ax,A.Eta,A.BalancedCorrelationScore,styles{m}, ...
        'LineWidth',2.0,'MarkerSize',7,'DisplayName',names{m});
end

submission_axes(ax);
paper_eta_axis(ax,eta_values);
submission_labels(ax,'Noise level, $\eta$','$\min(\bar{\rho}_C,\bar{\rho}_{\alpha})$');
ylim(ax,[-0.05 1.05]);
submission_legend(ax,'best');
save_figure_pair(fig,basepath);
close(fig);
end
