function [] = saveFigAndSvg(figfolder,figname)
% In figfolder, we assume there are two folders: MatlabFigs and SVGs. This
% function saves the current figure as figname in those folders in those
% formats.
%
% Inputs:
% - figfolder: string with the absolute or relative path to the folder of
% figures, which should contain two subfolders: "MatlabFigs" and "SVGs". In
% each of these, a .fig and .svg of the current figure will be saved,
% respectively. ASSUMES FIGFOLDER ENDS IN A \ OR /
% - figname: the name of the figure
%
% Adam Smoulder, 1/11/21

if ispc, slash = '\'; else, slash = '/'; end
saveas(gcf,[figfolder 'MatlabFigs' slash figname])
saveas(gcf,[figfolder 'SVGs' slash figname '.svg'])

end

