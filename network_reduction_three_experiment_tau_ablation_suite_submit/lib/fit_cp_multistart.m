function [best,run_table] = fit_cp_multistart(X,P,CP,seed0)
% Tensorlab standard CPD with CP rank P = R*q.
%
% Start 1 uses Tensorlab's high-level/default initialization.
% Starts 2,... explicitly use cpd_rnd.
% The selected run minimizes reconstruction NMSE on the OBSERVED tensor only.
% No clean truth is used to choose the CPD run.

    n_restarts = CP.n_restarts;

    rows = repmat(struct( ...
        'Restart',0, ...
        'Seed',0, ...
        'Initialization',"", ...
        'Success',false, ...
        'ObservedReconNMSE',Inf, ...
        'ElapsedSeconds',NaN, ...
        'ErrorMessage',""), ...
        n_restarts,1);

    results = cell(n_restarts,1);

    for restart = 1:n_restarts

        seed = seed0+restart;
        rng(seed,'twister');

        opts = struct;
        opts.Display = CP.display;

        % cpd() passes these to its optimization/refinement stages.
        opts.AlgorithmOptions.MaxIter = CP.max_iter;
        opts.AlgorithmOptions.TolFun = CP.tol_fun;
        opts.AlgorithmOptions.TolX = CP.tol_x;

        opts.RefinementOptions.MaxIter = CP.max_iter;
        opts.RefinementOptions.TolFun = CP.tol_fun;
        opts.RefinementOptions.TolX = CP.tol_x;

        if restart==1
            init_name = "Tensorlab default/high-level";
        else
            opts.Initialization = @cpd_rnd;
            init_name = "cpd_rnd";
        end

        rows(restart).Restart = restart;
        rows(restart).Seed = seed;
        rows(restart).Initialization = init_name;

        tic;

        try

            [Uhat,out] = ...
                cpd(X,P,opts);

            Xhat = ...
                real(cpdgen(Uhat));

            obs_nmse = ...
                sum((X(:)-Xhat(:)).^2) ...
                /max(sum(X(:).^2),eps);

            result = struct;
            result.Uhat = Uhat;
            result.output = out;
            result.Xhat = Xhat;
            result.observed_recon_nmse = obs_nmse;
            result.rank = P;

            results{restart} = result;

            rows(restart).Success = true;
            rows(restart).ObservedReconNMSE = obs_nmse;

        catch ME

            rows(restart).Success = false;
            rows(restart).ObservedReconNMSE = Inf;
            rows(restart).ErrorMessage = string(ME.message);
        end

        rows(restart).ElapsedSeconds = toc;
    end

    run_table = ...
        struct2table(rows);

    eligible = ...
        find( ...
        run_table.Success ...
        &isfinite(run_table.ObservedReconNMSE));

    if isempty(eligible)

        disp(run_table);

        error([ ...
            'All Tensorlab CPD starts failed. ' ...
            'Check cfg.tensorlabPath and Tensorlab CPD compatibility.']);
    end

    [~,k] = ...
        min(run_table.ObservedReconNMSE(eligible));

    ibest = eligible(k);

    best = results{ibest};
end
