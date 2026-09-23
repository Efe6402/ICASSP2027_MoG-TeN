function out = factor_correlations_balanced( ...
    Ctrue,atrue,Cest,aest,allow_sign_flip)
%FACTOR_CORRELATIONS_BALANCED
% Match estimated factors directly to the CLEAN factors by maximizing a
% balanced correlation criterion.
%
% For a candidate component permutation:
%
%   rho_C     = mean_r corr(vec(C_r^true), vec(C_r^est))
%   rho_alpha = mean_r corr(alpha_r^true, alpha_r^est)
%
% where corr is centered Pearson correlation.
%
% The BALANCED factor score is
%
%   score = min(rho_C,rho_alpha).
%
% This is intentionally conservative: one excellent factor family cannot
% compensate for a poor one.  The permutation maximizing this score is used.
%
% For BTD only, allow_sign_flip=true resolves the simultaneous
% (C_r,alpha_r) -> (-C_r,-alpha_r) ambiguity.

R = size(Ctrue,3);

assert(size(Cest,3)==R, ...
    'Cest must have the same number of components as Ctrue.');
assert(size(atrue,2)==R && size(aest,2)==R, ...
    'Temporal factor matrices must both have R columns.');

if R<=8
    P = perms(1:R);
else
    error(['Exact correlation-based matching currently supports R<=8. ' ...
           'These experiments use R=3.']);
end

best_score = -Inf;
best_average = -Inf;
best_worst = -Inf;
best_perm = 1:R;
best_signs = ones(R,1);
best_C = nan(R,1);
best_A = nan(R,1);

for pidx = 1:size(P,1)

    p = P(pidx,:);

    cvals = nan(R,1);
    avals = nan(R,1);
    signs = ones(R,1);

    for r = 1:R

        Ce = Cest(:,:,p(r));
        ae = aest(:,p(r));

        c0 = centered_corr_local(Ctrue(:,:,r),Ce);
        a0 = centered_corr_local(atrue(:,r),ae);

        if allow_sign_flip && isfinite(c0) && isfinite(a0)

            score_plus = min(c0,a0);
            score_minus = min(-c0,-a0);

            if score_minus > score_plus
                c0 = -c0;
                a0 = -a0;
                signs(r) = -1;
            end
        end

        cvals(r) = c0;
        avals(r) = a0;
    end

    if ~all(isfinite(cvals)) || ~all(isfinite(avals))
        continue;
    end

    rhoC = mean(cvals);
    rhoA = mean(avals);

    score = min(rhoC,rhoA);
    average_score = 0.5*(rhoC+rhoA);
    worst_individual = min([cvals;avals]);

    improve = false;

    if score > best_score + 1e-12
        improve = true;
    elseif abs(score-best_score)<=1e-12
        if average_score > best_average + 1e-12
            improve = true;
        elseif abs(average_score-best_average)<=1e-12 ...
                && worst_individual > best_worst + 1e-12
            improve = true;
        end
    end

    if improve
        best_score = score;
        best_average = average_score;
        best_worst = worst_individual;
        best_perm = p;
        best_signs = signs;
        best_C = cvals;
        best_A = avals;
    end
end

if ~isfinite(best_score)
    error('Could not compute a finite clean-factor correlation score.');
end

Caligned = zeros(size(Ctrue));
aaligned = zeros(size(atrue));

for r = 1:R

    Ce = Cest(:,:,best_perm(r));
    ae = aest(:,best_perm(r));

    if best_signs(r)<0
        Ce = -Ce;
        ae = -ae;
    end

    Caligned(:,:,r) = Ce;
    aaligned(:,r) = ae;
end

out.C_each = best_C;
out.alpha_each = best_A;

out.C_mean = mean(best_C);
out.alpha_mean = mean(best_A);

out.balanced_score = best_score;
out.average_score = best_average;
out.worst_individual_correlation = best_worst;

out.C_aligned = Caligned;
out.alpha_aligned = aaligned;
out.permutation = best_perm;
out.signs = best_signs;
end


function c = centered_corr_local(x,y)

x = real(x(:));
y = real(y(:));

x = x-mean(x);
y = y-mean(y);

nx = norm(x);
ny = norm(y);

if nx<=eps || ny<=eps
    c = NaN;
else
    c = (x'*y)/(nx*ny);
end

c = max(-1,min(1,real(c)));
end
