function configs = build_joint_hyperparameter_design(SEARCH,zeta)
% Positive joint search over gamma, tau, mu, lambda.
% tau is strictly positive in every configuration.

    anchors = [
        1e-4    1e-3    1e-2    1e-4
        3e-4    3e-3    3e-2    3e-4
        1e-3    1e-2    1e-1    1e-3
        3e-3    3e-2    3e-1    3e-3
        1e-2    1e-1    1       1e-2
        3e-2    2e-1    3       3e-2
        1e-1    5e-1    10      1e-1
        2e-1    1       20      3e-1
        3e-1    1.5     30      5e-1
        5e-1    2       50      1
    ];

    rng(SEARCH.seed_design,'twister');

    nr = SEARCH.n_random_design;

    ranges = [
        SEARCH.log10_gamma
        SEARCH.log10_tau
        SEARCH.log10_mu
        SEARCH.log10_lambda
    ];

    U = ...
        stratified_unit_design(nr,4);

    random_part = zeros(nr,4);

    for j = 1:4

        lo = ranges(j,1);
        hi = ranges(j,2);

        random_part(:,j) = ...
            10.^(lo+(hi-lo)*U(:,j));
    end

    P = [anchors;random_part];

    % Avoid effectively-zero overlap penalties in this dedicated lambda ablation.
    P = P(P(:,4)>=1e-3,:);

    assert(all(P(:)>0), ...
        'All regularization parameters must be strictly positive.');

    [~,ia] = unique(P,'rows','stable');
    P = P(sort(ia),:);

    n = size(P,1);

    configs = table( ...
        (1:n).', ...
        P(:,1),P(:,2),P(:,3),P(:,4), ...
        repmat(zeta,n,1), ...
        'VariableNames', { ...
        'ConfigID','gamma','tau','mu','lambda','zeta'});
end


function U = ...
    stratified_unit_design(n,d)

    U = zeros(n,d);

    for j = 1:d

        base = ...
            ((0:n-1)' ...
            +rand(n,1))/n;

        U(:,j) = ...
            base(randperm(n));
    end
end




%% ========================================================================
% PROPOSED SOLVER (NUCLEAR TERM ACTIVE THROUGH hyp.tau > 0)
% =========================================================================
