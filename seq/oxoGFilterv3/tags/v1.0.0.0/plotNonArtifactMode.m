function [] = plotNonArtifactMode(outputDir, uname, M_nA, acsY)
% TODO: Documentation
figure;
[N,xb,yb,h] = hist2d(M_nA.foxog, M_nA.alt_read_count,0:.01:1, 1:max(acsY));
hold;

% hh is the hist2d axes
hh = get(h(1), 'Parent');

% handle to the left histogram
hl = get(h(3), 'Parent');

text(.05, .9, [{'No C>A or G>T'};{['n = ' num2str(length(M_nA.foxog))]};], 'FontSize', 14, 'Units', 'Normalized', 'Parent', hh)
title(hh, [uname], 'FontSize', 16)
set(gca, 'FontSize', 14)
ylabel(hl, 'alt count', 'FontSize', 16)
xlabel('F_{oxoG}', 'FontSize', 16)
line([0 1], [0 0], 'Color', 'k')
print('-dpng', [outputDir uname '_NonArtifactMode_hist2d.png']);

close(gcf)
