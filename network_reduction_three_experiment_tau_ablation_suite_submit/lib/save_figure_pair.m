function save_figure_pair(fig,basepath)
%SAVE_FIGURE_PAIR Save MATLAB, high-resolution raster, and vector versions.
%
% Outputs:
%   <basepath>.fig
%   <basepath>.png   (300 dpi)
%   <basepath>.pdf   (vector when exportgraphics is available)

savefig(fig,[basepath '.fig']);

try
    exportgraphics(fig,[basepath '.png'],'Resolution',300);
catch
    print(fig,[basepath '.png'],'-dpng','-r300');
end

try
    exportgraphics(fig,[basepath '.pdf'],'ContentType','vector');
catch
    try
        print(fig,[basepath '.pdf'],'-dpdf','-painters');
    catch
    end
end
end
