function [] = plotNonArtifactMode(outputDir, uname, M_nA, acsY)
% [] = plotNonArtifactMode(outputDir, uname, M_nA, acsY)
%   Produce orchestra plot (2D histogram) with non-artifact mode mutation 
%       information.
%   The plot is produced as a png with the name:
%       [outputDir uname '_NonArtifactMode_hist2d.png']
%
% outputDir -- location to place generated png file.
%   
% uname -- case name (or etxt string to be used in title and png files).
%
% M_nA -- mutations, assumed to be non-artifact mode.
%
% acsY -- allele counts to be plotted.
%

figure;
[N,xb,yb,h] = hist2d(M_nA.foxog, M_nA.alt_read_count,0:.02:1, 1:max(acsY));
hold;
text(.05, .9, [{'Not OxoG Mode'};{['n = ' num2str(length(M_nA.foxog))]};], 'FontSize', 14, 'Units', 'Normalized', 'Parent', h(1))
title(h(1), [uname], 'FontSize', 14, 'Interpreter', 'None')
set(gca, 'FontSize', 14)
ylabel(h(3), 'alt count', 'FontSize', 16)
xlabel('F_{oxoG}', 'FontSize', 16,'Interpreter','tex')
saveas(gcf, [outputDir uname '_NonArtifactMode_hist2d.png']);
close(gcf)
