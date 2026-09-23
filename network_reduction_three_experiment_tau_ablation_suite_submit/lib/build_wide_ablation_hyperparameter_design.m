function configs = build_wide_ablation_hyperparameter_design(SEARCH,zeta)
%BUILD_WIDE_ABLATION_HYPERPARAMETER_DESIGN
% Broader strictly-positive search for the nuclear model.
%
% Columns are [gamma, tau, mu, lambda].  tau is NEVER zero here; tau=0 is
% created only in the paired ablation after a nuclear model has been selected.

anchors = [
    1e-5    1e-3    1e-3    1e-5
    3e-5    3e-3    3e-3    3e-5
    1e-4    1e-2    1e-2    1e-4
    3e-4    3e-2    3e-2    3e-4
    1e-3    1e-2    1e-1    1e-3
    3e-3    3e-2    3e-1    3e-3
    1e-2    1e-1    1       1e-2
    3e-2    2e-1    3       3e-2
    1e-1    5e-1    10      1e-1
    2e-1    1       30      3e-1
    5e-1    2       60      8e-1
    1       3       100     2
];

rng(SEARCH.seed_design,'twister');

nr = SEARCH.n_random_design;

ranges = [
    SEARCH.log10_gamma
    SEARCH.log10_tau
    SEARCH.log10_mu
    SEARCH.log10_lambda
];

U = stratified_unit_design_local(nr,4);
random_part = zeros(nr,4);

for j = 1:4
    lo = ranges(j,1);
    hi = ranges(j,2);
    random_part(:,j) = 10.^(lo+(hi-lo)*U(:,j));
end

P = [anchors; random_part];

% Keep the nuclear branch genuinely active.
P = P(P(:,2) >= SEARCH.min_tau_positive,:);

assert(all(P(:)>0), ...
    'All full-model regularization parameters must be strictly positive.');

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


function U = stratified_unit_design_local(n,d)

U = zeros(n,d);

for j = 1:d
    base = ((0:n-1)' + rand(n,1))/n;
    U(:,j) = base(randperm(n));
end
end
