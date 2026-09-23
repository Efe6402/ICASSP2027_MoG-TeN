function lgd = submission_legend(ax,location)
%SUBMISSION_LEGEND Compact academic legend used in submission figures.
%
% ------------------------- STYLE TUNING -------------------------
LEGEND_SIZE = 11; % <-- legend font size
% ---------------------------------------------------------------

if nargin<2 || isempty(location)
    location = 'best';
end

lgd = legend(ax, ...
    'Location',location, ...
    'Interpreter','latex', ...
    'FontSize',LEGEND_SIZE, ...
    'Box','off');
end
