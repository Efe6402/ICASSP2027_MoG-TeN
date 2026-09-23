function [X,Ccorrupt,NoiseTensor,actual_eta,actual_snr_db] = ...
    add_template_noise_at_eta(Cclean_raw,alpha_raw,Hbase,eta)
% EXACT SAME template-noise addition as the supplied script.
%
% Add noise DIRECTLY TO THE TEMPLATES while preserving unsigned graphs.
% First build a generic per-template direction proportional to ||C_r||_F.
% Then use ONE common scale so the induced tensor perturbation satisfies
%
%   ||X-Xclean||_F / ||Xclean||_F = eta.

    [N,~,R] = size(Cclean_raw);

    Delta0 = zeros(size(Cclean_raw));

    for r = 1:R

        Hr = Hbase(:,:,r);

        Hr = ...
            max((Hr+Hr')/2,0);

        Hr(1:N+1:end) = 0;

        Delta0(:,:,r) = ...
            norm(Cclean_raw(:,:,r),'fro') ...
            *Hr/max(norm(Hr,'fro'),eps);
    end

    Xclean = ...
        reconstruct_tensor(Cclean_raw,alpha_raw);

    Noise0 = ...
        reconstruct_tensor(Delta0,alpha_raw);

    scale = ...
        eta*norm(Xclean(:)) ...
        /max(norm(Noise0(:)),eps);

    Delta = scale*Delta0;

    Ccorrupt = Cclean_raw+Delta;

    % Numerical enforcement of unsigned symmetric hollow templates.
    for r = 1:R

        Cr = ...
            max((Ccorrupt(:,:,r)+Ccorrupt(:,:,r)')/2,0);

        Cr(1:N+1:end) = 0;

        Ccorrupt(:,:,r) = Cr;
    end

    X = ...
        reconstruct_tensor(Ccorrupt,alpha_raw);

    NoiseTensor = X-Xclean;

    actual_eta = ...
        norm(NoiseTensor(:)) ...
        /max(norm(Xclean(:)),eps);

    if norm(NoiseTensor(:))<=eps
        actual_snr_db = Inf;
    else
        actual_snr_db = ...
            20*log10( ...
            norm(Xclean(:))/norm(NoiseTensor(:)));
    end
end
