function [] = plotArtifactMode(outputDir, uname, M_A, fdrStruct, foxogX, acsY)
% TODO: Documentation
%
% Patch has been altered to produce a line.  Mutations to the right of the
%  line have been filtered.
%
% Putting [] for foxogX will cause this script to skip the
%  patch plotting in text blurb.

isPlottingPatch = 1;
if isempty(foxogX) || isempty(fdrStruct)
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

N_mut = length(fdrStruct.cut);
filteredMutationsString = ['Filtered Mutations: ' num2str(sum(fdrStruct.cut))];
filteredMutationsStringLine2 = ['  (' num2str(100*sum(fdrStruct.cut)/size(M_A.foxog,1) ,'%.1f')  '% of OxoG Mode)'];
filteredMutationsStringLine3 = ['  (' num2str(100*sum(fdrStruct.cut)/N_mut ,'%.1f') '% of all)'];

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
