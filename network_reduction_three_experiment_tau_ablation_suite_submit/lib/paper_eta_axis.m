function paper_eta_axis(ax,eta_values)
%PAPER_ETA_AXIS Clear logarithmic eta axis for geometric noise sweeps.
%
% The eta grid used in these experiments is approximately geometric.
% A log-x axis therefore both represents the sweep naturally and prevents
% the small-eta labels from colliding.

if nargin<1 || isempty(ax)
    ax = gca;
end

eta_values = eta_values(:).';

set(ax,'XScale','log');

n = numel(eta_values);

if n>=8
    idx = [1 3 5 7 n];
elseif n>=5
    idx = unique(round(linspace(1,n,5)));
else
    idx = 1:n;
end

ticks = eta_values(idx);

xticks(ax,ticks);
xticklabels(ax,arrayfun(@(x)sprintf('%.3g',x),ticks,'UniformOutput',false));

xlim(ax,[eta_values(1)/1.08, eta_values(end)*1.08]);
end
