function metrics = evaluate_tensor_reconstruction_only(X,Xclean,Xhat)
% Tensor-only metrics for standard CPD.
%
% Standard CPD with P=R*q returns P independent rank-one terms rather than
% R block components.  We intentionally do not cluster/group those terms
% into artificial C_r/alpha_r pairs.  Component-level fields are therefore
% NaN so the summary table remains schema-compatible with the other methods.

    metrics.noisy_recon_nmse = ...
        sum((X(:)-Xhat(:)).^2) ...
        /max(sum(X(:).^2),eps);

    metrics.clean_recon_nmse = ...
        sum((Xclean(:)-Xhat(:)).^2) ...
        /max(sum(Xclean(:).^2),eps);

    metrics.mean_C_relerr = NaN;
    metrics.mean_alpha_cosine = NaN;
    metrics.mean_tensor_component_relerr = NaN;
    metrics.min_alpha_cosine = NaN;
end
