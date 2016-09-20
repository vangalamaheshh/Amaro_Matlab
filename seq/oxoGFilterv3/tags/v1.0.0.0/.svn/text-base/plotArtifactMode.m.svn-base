function [] = plotArtifactMode(outputDir, uname, M_A, fdrStruct, foxogX, acsY)
% TODO: Documentation
% TODO: Summary histogram on the x-axis
figure;
[N,xb,yb,h] = hist2d(M_A.foxog , M_A.alt_read_count,0:.01:1, 1:max(acsY));

% hh is the hist2d axes
hh = get(h(1), 'Parent');

% handle to the left histogram
hl = get(h(3), 'Parent');

hold;
hp = patch([foxogX(1) foxogX foxogX(end) 1.1 1.1], [ 1 acsY max(acsY)+1 max(acsY)+1 1], [.8 .8 .8], 'Parent', hh);
set(hp, 'FaceAlpha', 0.3)
title(hh, [uname], 'FontSize', 16)
N_mut = length(fdrStruct.cut);
filteredMutationsString = ['Filtered Mutations: ' num2str(sum(fdrStruct.cut))];
filteredMutationsStringLine2 = ['  (' num2str(100*sum(fdrStruct.cut)/N_mut ,'%.1f') '/' num2str(100*sum(fdrStruct.cut)/size(M_A.foxog,1) ,'%.1f')  '% of all/OxoMode)'];
text(.05, .9, [{'C>A or G>T Only'};{['n = ' num2str(length(M_A.foxog))]};{filteredMutationsString};{filteredMutationsStringLine2}], 'FontSize', 14, 'Units', 'Normalized', 'Parent', hh)
set(hh, 'FontSize', 12)
ylabel(hl, 'alt count', 'FontSize', 16)
xlabel('F_{oxoG}', 'FontSize', 16)
line([0 1], [0 0], 'Color', 'k', 'Parent', hh)
print('-dpng', [outputDir uname '_ArtifactMode_hist2d_FDR.png']);

close(gcf)
