function lgd = paper_legend(ax,location)
%PAPER_LEGEND Consistent legend styling.

if nargin<1 || isempty(ax)
    ax = gca;
end

if nargin<2 || isempty(location)
    location = 'best';
end

lgd = legend(ax, ...
    'Location',location, ...
    'Interpreter','latex', ...
    'Box','off', ...
    'FontSize',11);
end
