function [best_result,best_metrics,best_corr,best_hyp, ...
    config_table,restart_table_all] = ...
    sweep_ours_hyperparameters_by_correlation( ...
    X,Xclean,Ctrue,atrue, ...
    R,L,V_L,ell_L, ...
    configs,SOLVER,n_restarts,seed0)
%SWEEP_OURS_HYPERPARAMETERS_BY_CORRELATION
%
% Within one hyperparameter tuple:
%   restart selection uses ONLY the ORIGINAL FITTED OBJECTIVE through
%   fit_ours_multistart().
%
% Across hyperparameter tuples:
%   CLEAN-X reconstruction error is NOT used for selection.
%
% Instead, after correlation-optimal component matching to CLEAN C_true and
% CLEAN alpha_true, select the tuple maximizing
%
%   min(mean C correlation, mean alpha correlation).
%
% Tie breaking:
%   1) maximum average of the two mean correlations,
%   2) maximum worst individual component/factor correlation,
%   3) minimum fitted objective.
%
% CLEAN-X NMSE is computed and saved only as a post-selection diagnostic.

nC = height(configs);

rows = repmat(struct( ...
    'ConfigID',0, ...
    'gamma',NaN, ...
    'tau',NaN, ...
    'mu',NaN, ...
    'lambda',NaN, ...
    'zeta',NaN, ...
    'MeanCCorrelation',-Inf, ...
    'MeanAlphaCorrelation',-Inf, ...
    'BalancedCorrelationScore',-Inf, ...
    'AverageCorrelationScore',-Inf, ...
    'WorstIndividualCorrelation',-Inf, ...
    'CleanReconNMSE',Inf, ...
    'ObservedReconNMSE',Inf, ...
    'MeanCRelErr',Inf, ...
    'MeanAlphaCosine',-Inf, ...
    'MeanComponentRelErr',Inf, ...
    'FinalObjective',Inf, ...
    'Converged',false, ...
    'Iterations',NaN), ...
    nC,1);

results = cell(nC,1);
metrics_all = cell(nC,1);
corr_all = cell(nC,1);

restart_table_all = table();

for cidx = 1:nC

    hyp.gamma = configs.gamma(cidx);
    hyp.tau = configs.tau(cidx);
    hyp.mu = configs.mu(cidx);
    hyp.lambda = configs.lambda(cidx);
    hyp.zeta = configs.zeta(cidx);

    % SAME restart seeds for every hyperparameter tuple.
    [candidate,candidate_runs] = fit_ours_multistart( ...
        X,R,L,V_L,ell_L,hyp,SOLVER,n_restarts,seed0);

    metrics = evaluate_solution_general( ...
        X,Xclean,candidate.C,candidate.alpha,Ctrue,atrue,false);

    corrdiag = factor_correlations_balanced( ...
        Ctrue,atrue,candidate.C,candidate.alpha,false);

    results{cidx} = candidate;
    metrics_all{cidx} = metrics;
    corr_all{cidx} = corrdiag;

    candidate_runs.ConfigID = ...
        repmat(configs.ConfigID(cidx),height(candidate_runs),1);
    candidate_runs.gamma = ...
        repmat(hyp.gamma,height(candidate_runs),1);
    candidate_runs.tau = ...
        repmat(hyp.tau,height(candidate_runs),1);
    candidate_runs.mu = ...
        repmat(hyp.mu,height(candidate_runs),1);
    candidate_runs.lambda = ...
        repmat(hyp.lambda,height(candidate_runs),1);

    restart_table_all = append_table( ...
        restart_table_all,candidate_runs);

    rows(cidx).ConfigID = configs.ConfigID(cidx);
    rows(cidx).gamma = hyp.gamma;
    rows(cidx).tau = hyp.tau;
    rows(cidx).mu = hyp.mu;
    rows(cidx).lambda = hyp.lambda;
    rows(cidx).zeta = hyp.zeta;

    rows(cidx).MeanCCorrelation = corrdiag.C_mean;
    rows(cidx).MeanAlphaCorrelation = corrdiag.alpha_mean;
    rows(cidx).BalancedCorrelationScore = corrdiag.balanced_score;
    rows(cidx).AverageCorrelationScore = corrdiag.average_score;
    rows(cidx).WorstIndividualCorrelation = ...
        corrdiag.worst_individual_correlation;

    rows(cidx).CleanReconNMSE = metrics.clean_recon_nmse;
    rows(cidx).ObservedReconNMSE = metrics.noisy_recon_nmse;
    rows(cidx).MeanCRelErr = metrics.mean_C_relerr;
    rows(cidx).MeanAlphaCosine = metrics.mean_alpha_cosine;
    rows(cidx).MeanComponentRelErr = ...
        metrics.mean_tensor_component_relerr;

    rows(cidx).FinalObjective = candidate.history.objective(end);
    rows(cidx).Converged = candidate.converged;
    rows(cidx).Iterations = candidate.iterations;
end

config_table = struct2table(rows);

eligible = find( ...
    isfinite(config_table.BalancedCorrelationScore) ...
    &isfinite(config_table.MeanCCorrelation) ...
    &isfinite(config_table.MeanAlphaCorrelation));

if isempty(eligible)
    error('All OURS hyperparameter configurations failed correlation scoring.');
end

% PRIMARY: maximize the worse of mean C and mean alpha correlations.
score = config_table.BalancedCorrelationScore(eligible);
best_score = max(score);
tol = 1e-12*max(1,abs(best_score));

tied = eligible(abs(score-best_score)<=tol);

% SECONDARY: maximize their average.
if numel(tied)>1
    avg = config_table.AverageCorrelationScore(tied);
    best_avg = max(avg);
    tol2 = 1e-12*max(1,abs(best_avg));
    tied = tied(abs(avg-best_avg)<=tol2);
end

% TERTIARY: maximize the worst individual matched correlation.
if numel(tied)>1
    w = config_table.WorstIndividualCorrelation(tied);
    best_w = max(w);
    tol3 = 1e-12*max(1,abs(best_w));
    tied = tied(abs(w-best_w)<=tol3);
end

% FINAL numerical tie-break: fitted objective, NOT clean-X.
if numel(tied)>1
    [~,k] = min(config_table.FinalObjective(tied));
    ibest = tied(k);
else
    ibest = tied(1);
end

best_result = results{ibest};
best_metrics = metrics_all{ibest};
best_corr = corr_all{ibest};

best_hyp.gamma = config_table.gamma(ibest);
best_hyp.tau = config_table.tau(ibest);
best_hyp.mu = config_table.mu(ibest);
best_hyp.lambda = config_table.lambda(ibest);
best_hyp.zeta = config_table.zeta(ibest);
best_hyp.ConfigID = config_table.ConfigID(ibest);

best_hyp.SelectionMetric = "balanced_clean_factor_correlation";
best_hyp.MeanCCorrelation = config_table.MeanCCorrelation(ibest);
best_hyp.MeanAlphaCorrelation = config_table.MeanAlphaCorrelation(ibest);
best_hyp.BalancedCorrelationScore = ...
    config_table.BalancedCorrelationScore(ibest);
end
