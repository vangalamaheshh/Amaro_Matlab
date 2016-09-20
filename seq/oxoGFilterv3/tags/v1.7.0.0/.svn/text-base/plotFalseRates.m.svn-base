function [] = plotFalseRates(pngFilename, titleToUse, ylabelToUse, Ns, falseVals)
%
% TODO: Documentation
%

figure;
xvals = [max(min(Ns)-10,0) max(Ns)+10];
line(logspace(0, log10(max(xvals)), 100), .01*logspace(0, log10(max(xvals)), 100), 'Color', 'g', 'LineWidth', 2)
hold on;
line(logspace(0, log10(max(xvals)), 100), .1*logspace(0, log10(max(xvals)), 100), 'Color', 'r', 'LineWidth', 2)
errorbar(Ns, falseVals(:,1),  falseVals(:,1)-falseVals(:,2), falseVals(:,3)-falseVals(:,1), 'o');
set(gca, 'XScale', 'log')
legend('One Percent', 'Ten Percent', 'Location', 'NorthWest');
title(titleToUse, 'Interpreter', 'None');
xlabel('Input Mutation Count');
ylabel(ylabelToUse);
print('-dpng', [pngFilename '.png']);

close(gcf)