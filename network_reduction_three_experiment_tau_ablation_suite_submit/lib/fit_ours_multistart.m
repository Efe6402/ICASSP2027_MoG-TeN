function [best,run_table] = fit_ours_multistart( ...
    X,R,L,V_L,ell_L,hyp,SOLVER,n_restarts,seed0)
% Multistart is strictly an optimization device.
% The selected run minimizes the ACTUAL fitted objective, never truth error.

    rows = repmat(struct( ...
        'Restart',0, ...
        'Seed',0, ...
        'Converged',false, ...
        'Iterations',0, ...
        'FinalObjective',Inf, ...
        'ObservedReconNMSE',Inf), ...
        n_restarts,1);

    results = cell(n_restarts,1);

    for restart = 1:n_restarts

        seed = seed0+restart;

        result = ...
            solve_network_reduction( ...
            X,R,L,V_L,ell_L,hyp,SOLVER,seed);

        results{restart} = result;

        rows(restart).Restart = restart;
        rows(restart).Seed = seed;
        rows(restart).Converged = result.converged;
        rows(restart).Iterations = result.iterations;
        rows(restart).FinalObjective = ...
            result.history.objective(end);
        rows(restart).ObservedReconNMSE = ...
            result.history.recon_nmse(end);
    end

    run_table = ...
        struct2table(rows);

    eligible = ...
        find( ...
        run_table.Converged ...
        &isfinite(run_table.FinalObjective));

    if isempty(eligible)
        eligible = ...
            find(isfinite(run_table.FinalObjective));
    end

    if isempty(eligible)
        error('All OURS multistart runs failed.');
    end

    [~,k] = ...
        min(run_table.FinalObjective(eligible));

    ibest = eligible(k);

    best = results{ibest};
end
