function C = generate_latent_group_templates( ...
    N,R,q,blocks,background_scale,target_fro,seed)

    rng(seed,'twister');

    C = zeros(N,N,R);

    for r = 1:R

        % ---------------------------------------------------------------
        % One-hot latent-group memberships.
        %
        % Every latent group is represented inside the dominant block and
        % also outside it.
        % ---------------------------------------------------------------
        labels = zeros(N,1);

        main_idx = blocks(r).indices;

        outside_idx = setdiff(1:N,main_idx,'stable');

        labels(main_idx) = ...
            balanced_random_labels( ...
            numel(main_idx), ...
            q);

        labels(outside_idx) = ...
            balanced_random_labels( ...
            numel(outside_idx), ...
            q);

        Z = zeros(N,q);

        for i = 1:N
            Z(i,labels(i)) = 1;
        end

        assert(rank(Z)==q, ...
            'Membership matrix Z is not full column rank.');

        % ---------------------------------------------------------------
        % Symmetric nonnegative hollow q-by-q interaction matrix.
        % Resample until numerically full rank.
        % ---------------------------------------------------------------
        W = make_full_rank_nonnegative_hollow_W(q);

        % ---------------------------------------------------------------
        % Node participation strengths.
        %
        % Main block: O(1)
        % Outside: faint but nonzero.
        %
        % Positive diagonal scaling preserves rank.
        % ---------------------------------------------------------------
        d = zeros(N,1);

        d(main_idx) = ...
            0.85 ...
            +0.30*rand(numel(main_idx),1);

        d(outside_idx) = ...
            background_scale ...
            *(0.70+0.60*rand(numel(outside_idx),1));

        D = diag(d);

        Cr = ...
            D ...
            *Z ...
            *W ...
            *Z' ...
            *D;

        Cr = (Cr+Cr')/2;

        % The diagonal is theoretically exactly zero because diag(W)=0.
        Cr(1:N+1:end) = 0;

        Cr = max(Cr,0);

        % Equal Frobenius energy across components/ranks.
        Cr = ...
            target_fro ...
            *Cr/max(norm(Cr,'fro'),eps);

        C(:,:,r) = Cr;

        nrank = numerical_rank(Cr);

        if nrank~=q
            error( ...
                'Option 1: expected rank %d but obtained %d.', ...
                q,nrank);
        end
    end
end

function W = make_full_rank_nonnegative_hollow_W(q)

    max_tries = 1000;

    for attempt = 1:max_tries

        A = ...
            0.25 ...
            +0.75*rand(q);

        A = triu(A,1);

        W = A+A';

        W(1:q+1:end) = 0;

        if numerical_rank(W)==q

            s = svd(W);

            if min(s)/max(s)>1e-4
                return;
            end
        end
    end

    error( ...
        'Could not generate a well-conditioned full-rank W.');
end


function labels = balanced_random_labels(n,q)

    labels = ...
        repmat((1:q).',ceil(n/q),1);

    labels = labels(1:n);

    labels = ...
        labels(randperm(n));
end


function r = numerical_rank(A)

    s = svd(A);

    if isempty(s) || max(s)==0
        r = 0;
        return;
    end

    tol = ...
        max(size(A)) ...
        *eps(max(s)) ...
        *10;

    r = sum(s>tol);
end

%% ========================================================================
% TRUTH CANONICALIZATION / BASIC OPERATORS
% =========================================================================
