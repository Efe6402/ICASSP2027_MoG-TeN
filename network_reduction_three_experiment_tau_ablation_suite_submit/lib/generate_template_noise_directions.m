function H = generate_template_noise_directions(N,R,seed)
% EXACT SAME template-noise direction construction as the supplied script.
%
% Generic high-rank UNSIGNED structural perturbations.
% Each H_r is symmetric, hollow, nonnegative, and unit Frobenius norm.

    rng(seed,'twister');

    H = zeros(N,N,R);

    for r = 1:R

        G = triu(abs(randn(N)),1);

        Hr = G+G';
        Hr(1:N+1:end) = 0;

        H(:,:,r) = ...
            Hr/max(norm(Hr,'fro'),eps);
    end
end
