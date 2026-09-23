%% RUNME_2_Q_SWEEP_WITH_TAU_ABLATION.M
%
% EXPERIMENT 2 -- TRUE-RANK SWEEP AT FIXED eta=0.316
% ================================================================
%
% q = [2 4 8 16]
%
% METHODS
%   Proposed
%   Proposed (tau=0), paired nuclear-norm ablation
%   BTD with L=[q q q]
%   CPD with P=3q
%
% OURS is independently RESWEPT at every q. Hyperparameters are selected by
%
%   max min(mean clean-C correlation, mean clean-alpha correlation).
%
% Then ONLY tau is set to zero for the ablation baseline.
%
% PRIMARY SUBMISSION FIGURE
%   submission_q_sweep_X_NMSE.png
%
% Additional factor-correlation diagnostics are saved separately.

% DONT FORGET TO REPLACE THE TENSORLAB PATH WITH THE USER'S OWN PATH 

clear; clc; close all;

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file), script_dir = pwd; else, script_dir = fileparts(this_file); end
addpath(genpath(fullfile(script_dir,'lib')));

CFG.tensorlabPath = ' '; %put your own tensorlab path 
assert(isfolder(CFG.tensorlabPath),'Tensorlab folder not found: %s',CFG.tensorlabPath);
addpath(genpath(CFG.tensorlabPath));
assert(exist('ll1','file')==2,'Tensorlab ll1.m not visible.');
assert(exist('cpd','file')==2,'Tensorlab cpd.m not visible.');

%% Data
DATA.N = 128;
DATA.T = 20;
DATA.R = 3;
DATA.q_values = [2 4 8 16];
DATA.target_C_fro = 10.0;
DATA.background_scale = 0.08;
DATA.seed_structure_base = 20262002;

TEMP.centers = [5.0 10.5 16.0];
TEMP.sigma = 3.0;

NOISE.eta = 0.3162277660;
NOISE.equivalent_snr_db = 10.0;
NOISE.seed_tuning = 20276001;
NOISE.seed_eval_base = 20277000;

EVAL.n_repetitions = 3;

%% OURS search
SEARCH.n_random_design = 10;
SEARCH.seed_design = 20260915;
SEARCH.log10_gamma  = [-4, log10(0.5)];
SEARCH.log10_tau    = [-3, log10(2.0)];
SEARCH.log10_mu     = [-2, log10(50)];
SEARCH.log10_lambda = [-3, 0];
SEARCH.min_lambda_positive = 1e-3;

TUNE.n_restarts = 3;
TUNE.seed_solver0 = 118501;

SOLVER.zeta = 2.0;
SOLVER.max_iter = 240;
SOLVER.min_iter = 15;
SOLVER.tol_primal = 1e-4;
SOLVER.tol_dual = 2e-4;
SOLVER.tol_obj = 1e-6;
SOLVER.alpha_tol = 1e-10;
SOLVER.range_tol = 1e-10;
SOLVER.eig_tol = 1e-12;
SOLVER.bisect_tol = 1e-12;
SOLVER.bisect_maxit = 110;
SOLVER.init_density = 0.15;
SOLVER.verbose_every = 0;

%% Evaluation fits
OURS_EVAL.n_restarts = 3;
OURS_EVAL.seed0 = 119001;

BTD.n_restarts = 2;
BTD.seed0 = 120001;
BTD.max_iter = 350;
BTD.tol_fun = 1e-10;
BTD.tol_x = 1e-8;
BTD.display = 0;

CP.n_restarts = 2;
CP.seed0 = 121001;
CP.max_iter = 350;
CP.tol_fun = 1e-10;
CP.tol_x = 1e-8;
CP.display = 0;

%% Output
timestamp = datestr(now,'yyyymmdd_HHMMSS');
root_dir = fullfile(script_dir,['results_02_q_tau_ablation_' timestamp]);
fig_dir = fullfile(root_dir,'submission');
diag_dir = fullfile(root_dir,'diagnostics');
submission_dir = fullfile(script_dir,'submission_figures');
mkdir(root_dir); mkdir(fig_dir); mkdir(diag_dir);
if ~isfolder(submission_dir), mkdir(submission_dir); end

try
    copyfile(this_file,fullfile(root_dir,'RUNME_2_q_sweep_with_tau_ablation.m'));
catch
end

%% Shared temporal truth / candidates / noise
alpha_true_raw = make_gaussian_temporal_factors( ...
    DATA.T,DATA.R,TEMP.centers,TEMP.sigma);

configs = build_qsweep_hyperparameter_design(SEARCH,SOLVER.zeta);
writetable(configs,fullfile(root_dir,'positive_tau_hyperparameter_candidates.csv'));

[L,V_L,ell_L] = temporal_laplacian(DATA.T);

H_tune = generate_template_noise_directions(DATA.N,DATA.R,NOISE.seed_tuning);
H_eval = cell(EVAL.n_repetitions,1);
for rep = 1:EVAL.n_repetitions
    H_eval{rep} = generate_template_noise_directions( ...
        DATA.N,DATA.R,NOISE.seed_eval_base+rep);
end

all_summary = table();
all_components = table();
all_hyp = table();

%% q sweep
for qidx = 1:numel(DATA.q_values)

    q = DATA.q_values(qidx);

    fprintf('\n============================================================\n');
    fprintf('Q SWEEP | q=%d\n',q);
    fprintf('============================================================\n');

    q_dir = fullfile(root_dir,sprintf('q_%02d',q));
    mkdir(q_dir);

    blocks = make_component_blocks(DATA.N,DATA.R);
    structure_seed = DATA.seed_structure_base + q;

    C_true_raw = generate_latent_group_templates( ...
        DATA.N,DATA.R,q,blocks, ...
        DATA.background_scale,DATA.target_C_fro,structure_seed);

    X_clean = reconstruct_tensor(C_true_raw,alpha_true_raw);

    [C_true,alpha_true,alpha_l2_before,canon_error] = ...
        canonicalize_truth(C_true_raw,alpha_true_raw);

    assert(canon_error<1e-10,'Truth canonicalization error at q=%d.',q);

    save(fullfile(q_dir,'clean_ground_truth.mat'), ...
        'q','C_true_raw','C_true','alpha_true_raw','alpha_true', ...
        'alpha_l2_before','canon_error','X_clean','-v7.3');

    %% Independent tuning observation
    [X_tune,C_corrupt_tune,Noise_tune,actual_eta_tune,actual_snr_tune] = ...
        add_template_noise_at_eta(C_true_raw,alpha_true_raw,H_tune,NOISE.eta);

    save(fullfile(q_dir,'tuning_observation.mat'), ...
        'X_tune','C_corrupt_tune','Noise_tune', ...
        'actual_eta_tune','actual_snr_tune','-v7.3');

    %% Positive-tau full model: resweep and select by component correlation
    seed_tune = TUNE.seed_solver0 + 1000*qidx;

    [ours_tuned,ours_tune_metrics,ours_tune_corr,ours_hyp, ...
        ours_config_table,ours_restart_table] = ...
        sweep_ours_hyperparameters_by_correlation( ...
        X_tune,X_clean,C_true,alpha_true, ...
        DATA.R,L,V_L,ell_L,configs,SOLVER,TUNE.n_restarts,seed_tune);

    writetable(ours_config_table,fullfile(q_dir,'positive_tau_config_search.csv'));
    writetable(ours_restart_table,fullfile(q_dir,'positive_tau_all_restarts.csv'));

    tau0_hyp = ours_hyp;
    tau0_hyp.tau = 0;

    all_hyp = [all_hyp; table( ...
        q,ours_hyp.ConfigID,ours_hyp.gamma,ours_hyp.tau,ours_hyp.mu, ...
        ours_hyp.lambda,ours_hyp.zeta,ours_tune_corr.C_mean, ...
        ours_tune_corr.alpha_mean,ours_tune_corr.balanced_score, ...
        'VariableNames',{'TrueRank','ConfigID','gamma','tau','mu','lambda','zeta', ...
        'TuningMeanCCorrelation','TuningMeanAlphaCorrelation', ...
        'TuningBalancedCorrelationScore'})]; %#ok<AGROW>

    q_summary = table();
    q_components = table();

    Lvec = q*ones(1,DATA.R);
    Pcp = DATA.R*q;

    %% Held-out evaluation
    for rep = 1:EVAL.n_repetitions

        [Xobs,C_corrupt_eval,Noise_eval,actual_eta,actual_snr_db] = ...
            add_template_noise_at_eta(C_true_raw,alpha_true_raw,H_eval{rep},NOISE.eta);

        seed_ours = OURS_EVAL.seed0 + 10000*qidx + 100*rep;

        [ours_result,ours_runs] = fit_ours_multistart( ...
            Xobs,DATA.R,L,V_L,ell_L,ours_hyp,SOLVER,OURS_EVAL.n_restarts,seed_ours);

        [tau0_result,tau0_runs] = fit_ours_multistart( ...
            Xobs,DATA.R,L,V_L,ell_L,tau0_hyp,SOLVER,OURS_EVAL.n_restarts,seed_ours);

        ours_metrics = evaluate_solution_general( ...
            Xobs,X_clean,ours_result.C,ours_result.alpha,C_true,alpha_true,false);
        ours_corr = factor_correlations_balanced( ...
            C_true,alpha_true,ours_result.C,ours_result.alpha,false);

        tau0_metrics = evaluate_solution_general( ...
            Xobs,X_clean,tau0_result.C,tau0_result.alpha,C_true,alpha_true,false);
        tau0_corr = factor_correlations_balanced( ...
            C_true,alpha_true,tau0_result.C,tau0_result.alpha,false);

        btd_result = fit_btd_multistart( ...
            Xobs,Lvec,BTD,BTD.seed0+10000*qidx+100*rep);
        btd_metrics = evaluate_solution_general( ...
            Xobs,X_clean,btd_result.C,btd_result.alpha,C_true,alpha_true,true);
        btd_corr = factor_correlations_balanced( ...
            C_true,alpha_true,btd_result.C,btd_result.alpha,true);

        cp_result = fit_cp_multistart( ...
            Xobs,Pcp,CP,CP.seed0+10000*qidx+100*rep);
        cp_metrics = evaluate_tensor_reconstruction_only(Xobs,X_clean,cp_result.Xhat);

        q_summary = append_summary(q_summary,q,rep,actual_eta,'proposed', ...
            ours_metrics.clean_recon_nmse,ours_metrics.noisy_recon_nmse, ...
            ours_corr.C_mean,ours_corr.alpha_mean,ours_corr.balanced_score);

        q_summary = append_summary(q_summary,q,rep,actual_eta,'tau0', ...
            tau0_metrics.clean_recon_nmse,tau0_metrics.noisy_recon_nmse, ...
            tau0_corr.C_mean,tau0_corr.alpha_mean,tau0_corr.balanced_score);

        q_summary = append_summary(q_summary,q,rep,actual_eta,'btd', ...
            btd_metrics.clean_recon_nmse,btd_metrics.noisy_recon_nmse, ...
            btd_corr.C_mean,btd_corr.alpha_mean,btd_corr.balanced_score);

        q_summary = append_summary(q_summary,q,rep,actual_eta,'cp', ...
            cp_metrics.clean_recon_nmse,cp_metrics.noisy_recon_nmse,NaN,NaN,NaN);

        for r = 1:DATA.R
            q_components = append_component(q_components,q,rep,'proposed',r, ...
                ours_corr.C_each(r),ours_corr.alpha_each(r));
            q_components = append_component(q_components,q,rep,'tau0',r, ...
                tau0_corr.C_each(r),tau0_corr.alpha_each(r));
            q_components = append_component(q_components,q,rep,'btd',r, ...
                btd_corr.C_each(r),btd_corr.alpha_each(r));
        end

        save(fullfile(q_dir,sprintf('heldout_rep_%02d.mat',rep)), ...
            'q','rep','X_clean','Xobs','C_corrupt_eval','Noise_eval', ...
            'ours_result','ours_runs','ours_metrics','ours_corr','ours_hyp', ...
            'tau0_result','tau0_runs','tau0_metrics','tau0_corr','tau0_hyp', ...
            'btd_result','btd_metrics','btd_corr','cp_result','cp_metrics','-v7.3');
    end

    writetable(q_summary,fullfile(q_dir,'heldout_summary.csv'));
    writetable(q_components,fullfile(q_dir,'heldout_component_correlations.csv'));

    all_summary = [all_summary;q_summary]; %#ok<AGROW>
    all_components = [all_components;q_components]; %#ok<AGROW>

    writetable(all_summary,fullfile(root_dir,'all_q_summary_partial.csv'));
    writetable(all_components,fullfile(root_dir,'all_q_component_correlations_partial.csv'));
    writetable(all_hyp,fullfile(root_dir,'selected_hyperparameters_partial.csv'));
end

%% Aggregate
writetable(all_summary,fullfile(root_dir,'all_q_heldout_summary.csv'));
writetable(all_components,fullfile(root_dir,'all_q_component_correlations.csv'));
writetable(all_hyp,fullfile(root_dir,'selected_hyperparameters_by_q.csv'));

S = aggregate_by_q(all_summary,DATA.q_values);
writetable(S,fullfile(root_dir,'all_q_summary_mean_std.csv'));

%% Submission figure
fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

plot_q_method(ax,S,"proposed",'Proposed','-o');
plot_q_method(ax,S,"tau0",'Proposed ($\tau=0$)','--s');
plot_q_method(ax,S,"btd",'BTD','-.^');
plot_q_method(ax,S,"cp",'CPD',':d');

submission_axes(ax);
xticks(ax,DATA.q_values);
xticklabels(ax,string(DATA.q_values));
xlim(ax,[1.5 16.5]);
submission_labels(ax,'True template rank, $q$','Tensor $\mathcal{X}$ Reconstruction NMSE');
submission_legend(ax,'best');

% Intentionally no title.
save_figure_pair(fig,fullfile(fig_dir,'submission_q_sweep_X_NMSE'));
save_figure_pair(fig,fullfile(submission_dir,'submission_q_sweep_X_NMSE'));
close(fig);

%% Diagnostics
plot_q_corr(S,DATA.q_values,'C',fullfile(diag_dir,'diagnostic_C_correlation_vs_q'));
plot_q_corr(S,DATA.q_values,'alpha',fullfile(diag_dir,'diagnostic_alpha_correlation_vs_q'));
plot_q_balanced(S,DATA.q_values,fullfile(diag_dir,'diagnostic_balanced_factor_score_vs_q'));

fprintf('\nTASK 2 COMPLETE\nResults: %s\n',root_dir);


%% ========================================================================
% LOCAL HELPERS
% =========================================================================

function T = append_summary(T,q,rep,eta,method,clean_nmse,observed_nmse,ccorr,acorr,score)
row = table(q,rep,eta,string(method),clean_nmse,observed_nmse,ccorr,acorr,score, ...
    'VariableNames',{'TrueRank','Repetition','Eta','Method','CleanReconNMSE', ...
    'ObservedReconNMSE','MeanCCorrelation','MeanAlphaCorrelation', ...
    'BalancedCorrelationScore'});
T = [T;row];
end


function T = append_component(T,q,rep,method,r,ccorr,acorr)
row = table(q,rep,string(method),r,ccorr,acorr, ...
    'VariableNames',{'TrueRank','Repetition','Method','Component', ...
    'CCorrelation','AlphaCorrelation'});
T = [T;row];
end


function S = aggregate_by_q(T,q_values)
methods = unique(T.Method,'stable');
S = table();

for m = 1:numel(methods)
    for q = q_values
        A = T(T.Method==methods(m) & T.TrueRank==q,:);
        row = table(methods(m),q, ...
            mean(A.CleanReconNMSE,'omitnan'),std(A.CleanReconNMSE,'omitnan'), ...
            mean(A.MeanCCorrelation,'omitnan'),std(A.MeanCCorrelation,'omitnan'), ...
            mean(A.MeanAlphaCorrelation,'omitnan'),std(A.MeanAlphaCorrelation,'omitnan'), ...
            mean(A.BalancedCorrelationScore,'omitnan'),std(A.BalancedCorrelationScore,'omitnan'), ...
            'VariableNames',{'Method','TrueRank','CleanReconNMSE_Mean', ...
            'CleanReconNMSE_Std','CCorrelation_Mean','CCorrelation_Std', ...
            'AlphaCorrelation_Mean','AlphaCorrelation_Std', ...
            'BalancedScore_Mean','BalancedScore_Std'});
        S = [S;row];
    end
end
end


function plot_q_method(ax,S,method,label,style)
A = sortrows(S(S.Method==string(method),:),'TrueRank');
errorbar(ax,A.TrueRank,A.CleanReconNMSE_Mean,A.CleanReconNMSE_Std,style, ...
    'LineWidth',2.2,'MarkerSize',7,'CapSize',7,'DisplayName',label);
end


function plot_q_corr(S,q_values,which_factor,basepath)

methods = ["proposed","tau0","btd"];
names = {'Proposed','Proposed ($\tau=0$)','BTD'};
styles = {'-o','--s','-.^'};

fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

for m = 1:numel(methods)
    A = sortrows(S(S.Method==methods(m),:),'TrueRank');

    if strcmp(which_factor,'C')
        yy = A.CCorrelation_Mean; ee = A.CCorrelation_Std;
    else
        yy = A.AlphaCorrelation_Mean; ee = A.AlphaCorrelation_Std;
    end

    errorbar(ax,A.TrueRank,yy,ee,styles{m}, ...
        'LineWidth',1.8,'MarkerSize',7,'CapSize',7,'DisplayName',names{m});
end

submission_axes(ax);
xticks(ax,q_values); xticklabels(ax,string(q_values));

if strcmp(which_factor,'C')
    submission_labels(ax,'True template rank, $q$','Mean $C_r$ Correlation');
else
    submission_labels(ax,'True template rank, $q$','Mean $\alpha_r$ Correlation');
end

ylim(ax,[-0.05 1.05]);
submission_legend(ax,'best');
save_figure_pair(fig,basepath);
close(fig);
end


function plot_q_balanced(S,q_values,basepath)

methods = ["proposed","tau0","btd"];
names = {'Proposed','Proposed ($\tau=0$)','BTD'};
styles = {'-o','--s','-.^'};

fig = figure('Visible','off','Color','w','Position',[100 100 900 590]);
ax = axes(fig); hold(ax,'on');

for m = 1:numel(methods)
    A = sortrows(S(S.Method==methods(m),:),'TrueRank');
    errorbar(ax,A.TrueRank,A.BalancedScore_Mean,A.BalancedScore_Std,styles{m}, ...
        'LineWidth',1.8,'MarkerSize',7,'CapSize',7,'DisplayName',names{m});
end

submission_axes(ax);
xticks(ax,q_values); xticklabels(ax,string(q_values));
submission_labels(ax,'True template rank, $q$','$\min(\bar{\rho}_C,\bar{\rho}_{\alpha})$');
ylim(ax,[-0.05 1.05]);
submission_legend(ax,'best');
save_figure_pair(fig,basepath);
close(fig);
end
