function [] = plotNonArtifactMode(outputDir, uname, M_nA, acsY)
% TODO: Documentation
figure;
[N,xb,yb,h] = hist2d(M_nA.foxog, M_nA.alt_read_count,0:.01:1, 1:max(acsY));
hold;

text(.05, .9, [{'No C>A or G>T'};{['n = ' num2str(length(M_nA.foxog))]};], 'FontSize', 14, 'Units', 'Normalized', 'Parent', h(1))
title(h(1), [uname], 'FontSize', 16)
set(gca, 'FontSize', 14)
ylabel(h(3), 'alt count', 'FontSize', 16)
xlabel('F_{oxoG}', 'FontSize', 16)
print('-dpng', [outputDir uname '_NonArtifactMode_hist2d.png']);

close(gcf)
