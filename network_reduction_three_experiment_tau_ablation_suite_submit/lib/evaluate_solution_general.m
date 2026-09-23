function metrics = evaluate_solution_general( ...
    X,Xclean,Cest,aest,Ctrue,atrue,allow_sign_flip)
% Match complete tensor blocks first, then score clean C and alpha.
%
% For BTD only:
%   (C,a) and (-C,-a) are equivalent and sign is resolved after matching.

    Xhat = ...
        reconstruct_tensor(Cest,aest);

    metrics.noisy_recon_nmse = ...
        sum((X(:)-Xhat(:)).^2) ...
        /max(sum(X(:).^2),eps);

    metrics.clean_recon_nmse = ...
        sum((Xclean(:)-Xhat(:)).^2) ...
        /max(sum(Xclean(:).^2),eps);

    R = size(Ctrue,3);

    cost = zeros(R,R);

    for rt = 1:R

        Gtrue = ...
            component_tensor( ...
            Ctrue(:,:,rt),atrue(:,rt));

        normG = ...
            max(norm(Gtrue(:)),eps);

        for re = 1:R

            Gest = ...
                component_tensor( ...
                Cest(:,:,re),aest(:,re));

            cost(rt,re) = ...
                norm(Gest(:)-Gtrue(:)) ...
                /normG;
        end
    end

    perm = ...
        best_permutation(cost);

    matched_comp = zeros(R,1);
    matched_C = zeros(R,1);
    matched_alpha_cos = zeros(R,1);
    signs = ones(R,1);

    for rt = 1:R

        re = perm(rt);

        Ce = Cest(:,:,re);
        ae = aest(:,re);

        if allow_sign_flip

            cc = ...
                safe_cosine( ...
                atrue(:,rt),ae);

            if isfinite(cc) && cc<0

                Ce = -Ce;
                ae = -ae;

                signs(rt) = -1;
            end
        end

        Gtrue = ...
            component_tensor( ...
            Ctrue(:,:,rt),atrue(:,rt));

        Gest = ...
            component_tensor(Ce,ae);

        matched_comp(rt) = ...
            norm(Gest(:)-Gtrue(:)) ...
            /max(norm(Gtrue(:)),eps);

        matched_C(rt) = ...
            norm(Ce-Ctrue(:,:,rt),'fro') ...
            /max(norm(Ctrue(:,:,rt),'fro'),eps);

        matched_alpha_cos(rt) = ...
            safe_cosine( ...
            atrue(:,rt),ae);
    end

    metrics.permutation = perm;
    metrics.signs = signs;
    metrics.matching_cost_matrix = cost;

    metrics.component_relerr = matched_comp;
    metrics.C_relerr = matched_C;
    metrics.alpha_cosine = matched_alpha_cos;

    metrics.mean_tensor_component_relerr = ...
        mean(matched_comp);

    metrics.mean_C_relerr = ...
        mean(matched_C);

    metrics.mean_alpha_cosine = ...
        mean(matched_alpha_cos);

    metrics.min_alpha_cosine = ...
        min(matched_alpha_cos);
end


function G = component_tensor(C,a)

    N = size(C,1);
    T = length(a);

    G = zeros(N,N,T);

    for t = 1:T
        G(:,:,t) = a(t)*C;
    end
end


function perm_best = ...
    best_permutation(cost)

    R = size(cost,1);

    if R<=8

        P = perms(1:R);

        best = Inf;
        perm_best = 1:R;

        for i = 1:size(P,1)

            p = P(i,:);

            s = 0;

            for r = 1:R
                s = ...
                    s+cost(r,p(r));
            end

            if s<best
                best = s;
                perm_best = p;
            end
        end

    else

        perm_best = zeros(1,R);

        available = true(1,R);

        for r = 1:R

            vals = cost(r,:);

            vals(~available) = Inf;

            [~,j] = min(vals);

            perm_best(r) = j;

            available(j) = false;
        end
    end
end


function c = safe_cosine(x,y)

    nx = norm(x(:));
    ny = norm(y(:));

    if nx<=eps || ny<=eps
        c = NaN;
    else
        c = ...
            (x(:)'*y(:))/(nx*ny);
    end
end


%% ========================================================================
% CLEAN DATA GENERATORS AND PROPOSED SOLVER
% Copied from the supplied exact-rank comparison code.
% =========================================================================
