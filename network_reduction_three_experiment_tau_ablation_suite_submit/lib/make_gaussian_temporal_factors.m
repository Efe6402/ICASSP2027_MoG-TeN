function alpha = make_gaussian_temporal_factors(T,R,centers,sigma)
%MAKE_GAUSSIAN_TEMPORAL_FACTORS Gaussian temporal activations, column-sum normalized.
assert(numel(centers)==R,'centers must contain R entries.');
t = (1:T).';
alpha = zeros(T,R);
for r = 1:R
    alpha(:,r) = exp(-0.5*((t-centers(r))/sigma).^2);
    alpha(:,r) = alpha(:,r)/max(sum(alpha(:,r)),eps);
end
end
