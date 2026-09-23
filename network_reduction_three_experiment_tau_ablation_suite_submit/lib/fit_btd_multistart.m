function [best,run_table] = fit_btd_multistart(X,Lvec,BTD,seed0)
% Tensorlab LL1 decomposition with known L.
%
% Start 1 uses Tensorlab's default/principal initialization.
% Starts 2,... use ll1_rnd.
% Selection uses observed reconstruction NMSE only.

    n_restarts = BTD.n_restarts;

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

        opts.Display = BTD.display;
        opts.MaxIter = BTD.max_iter;
        opts.TolFun = BTD.tol_fun;
        opts.TolX = BTD.tol_x;
        opts.OutputFormat = 'btd';

        if restart==1
            init_name = "Tensorlab default/principal";
        else
            opts.Initialization = @ll1_rnd;
            init_name = "ll1_rnd";
        end

        rows(restart).Restart = restart;
        rows(restart).Seed = seed;
        rows(restart).Initialization = init_name;

        tic;

        try

            [Uhat,out] = ...
                ll1(X,Lvec,opts);

            Ubtd = ...
                ensure_ll1_btd_format(Uhat,Lvec);

            try
                Xhat = ll1gen(Ubtd);
            catch
                Xhat = ll1gen(Ubtd,Lvec);
            end

            Xhat = real(Xhat);

            [Cest,aest] = ...
                extract_ll1_btd_components(Ubtd);

            obs_nmse = ...
                sum((X(:)-Xhat(:)).^2) ...
                /max(sum(X(:).^2),eps);

            result = struct;

            result.Ubtd = Ubtd;
            result.output = out;
            result.C = Cest;
            result.alpha = aest;
            result.Xhat = Xhat;
            result.observed_recon_nmse = obs_nmse;

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
            'All Tensorlab LL1/BTD starts failed. ' ...
            'Check cfg.tensorlabPath and Tensorlab compatibility.']);
    end

    [~,k] = ...
        min(run_table.ObservedReconNMSE(eligible));

    ibest = eligible(k);

    best = results{ibest};
end


function Ubtd = ensure_ll1_btd_format(U,Lvec)
% Tensorlab's default LL1 output is BTD format.
% This helper also accepts CPD-format output.

    if iscell(U) ...
            &&~isempty(U) ...
            &&iscell(U{1})

        Ubtd = U;

        return;
    end

    Ubtd = ...
        ll1convert(U,Lvec);
end


function [C,alpha] = extract_ll1_btd_components(Ubtd)
% Extract M_r=A_r*S_r*B_r' and temporal c_r.
% Canonicalize c_r to unit l2 norm, absorbing positive scale into M_r.

    R = numel(Ubtd);

    A1 = Ubtd{1}{1};
    c1 = Ubtd{1}{3};

    N = size(A1,1);
    T = numel(c1);

    C = zeros(N,N,R);
    alpha = zeros(T,R);

    for r = 1:R

        Ur = Ubtd{r};

        A = real(Ur{1});
        B = real(Ur{2});
        c = real(Ur{3}(:));
        S = real(Ur{4});

        if ndims(S)>2
            S = squeeze(S);
        end

        if isscalar(S)
            M = S*(A*B');
        else
            M = A*S*B';
        end

        sc = norm(c,2);

        if sc<=eps

            alpha(:,r) = c;
            C(:,:,r) = M;

        else

            alpha(:,r) = c/sc;
            C(:,:,r) = sc*M;
        end
    end
end
