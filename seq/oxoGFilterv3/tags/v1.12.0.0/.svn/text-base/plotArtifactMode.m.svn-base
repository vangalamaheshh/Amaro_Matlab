function [] = plotArtifactMode(outputDir, uname, M_A, cut, foxogX, acsY)
% [] = plotArtifactMode(outputDir, uname, M_A, cut, foxogX, acsY)
%   Produce orchestra plot (2D histogram) with artifact mode information 
%       and, optionally, a cut line.  All mutations in bins to right of 
%       the cutline are filtered.
%   The plot is produced as a png with the name:
%       [outputDir uname '_ArtifactMode_hist2d_FDR.png']
%
% outputDir -- location to place generated png file.
%   
% uname -- case name (or etxt string to be used in title and png files).
% 
% cut -- vector of true/false corresponding to the M_A structure array.
%
% foxogX and acsY -- x and y coordinates that define a cut line.  
%
% Mutations to the right of the line have been filtered.
%
% Putting [] for foxogX will cause this script to skip the
%  cut line plot.
% 
% See also plotNonArtifactMode

isPlottingPatch = 1;
if isempty(foxogX) || isempty(cut)
   disp(['No plotting a patch: ' uname]); 
   isPlottingPatch = 0;
end

figure;
[N,xb,yb,h] = hist2d(M_A.foxog , M_A.alt_read_count,0:.02:1, 1:max(acsY));

hold;

if isPlottingPatch
%     hp = patch([foxogX(1) foxogX foxogX(end) 1.1 1.1], [ 1 acsY max(acsY)+1 max(acsY)+1 1], [.8 .8 .8], 'Parent', h(1));
%     set(hp, 'FaceAlpha', 0.3)
    line([foxogX(1) foxogX ], [ 1 acsY ], 'Parent', h(1), 'LineWidth',1,'Color','k');
end
title(h(1), [uname], 'FontSize', 14, 'Interpreter', 'None')

N_mut = length(cut);
filteredMutationsString = ['Filtered Mutations: ' num2str(sum(cut))];
filteredMutationsStringLine2 = ['  (' num2str(100*sum(cut)/size(M_A.foxog,1) ,'%.1f')  '% of OxoG Mode)'];
filteredMutationsStringLine3 = ['  (' num2str(100*sum(cut)/N_mut ,'%.1f') '% of all)'];

text(.05, .85, [{'OxoG Mode'};{['n = ' num2str(length(M_A.foxog))]};{filteredMutationsString};{filteredMutationsStringLine2};{filteredMutationsStringLine3}], 'FontSize', 14, 'Units', 'Normalized', 'Parent', h(1))
set(h(1), 'FontSize', 12)
ylabel(h(3), 'alt count', 'FontSize', 16)
xlabel('F_{oxoG}', 'FontSize', 16)
line([0 1], [0 0], 'Color', 'k', 'Parent', h(1))
% print('-dpng', [outputDir uname '_ArtifactMode_hist2d_FDR.png']);
disp('Saving...')
saveas(gcf, [outputDir uname '_ArtifactMode_hist2d_FDR.png']);
disp('Done...')
close(gcf)
