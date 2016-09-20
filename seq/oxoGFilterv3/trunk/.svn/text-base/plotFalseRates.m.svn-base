function [] = plotFalseRates(pngFilename, titleToUse, ylabelToUse, Ns, falseVals)
%
% [] = plotFalseRates(pngFilename, titleToUse, ylabelToUse, Ns, falseVals)
%
%   Generates a plot with 95% confidence intervals for all cases given.
%
%       pngFilename -- name of the png to generate that contains the plot.
%
%       titleToUse -- The title on the plot.
%
%       ylabelToUse -- What should be on the ylabel
%
%       Ns -- nx1 matrix.  Number of mutations for each case.
%
%       falseVals -- nx3 matrix.  First column is the estimate.  Second 
%           column is the low CI.  Third column is the high CI.  Each row
%           is a case.
%

figure;
xvals = [max(min(Ns)-10,0) max(Ns)+10];
line(logspace(0, log10(max(xvals)), 100), .01*logspace(0, log10(max(xvals)), 100), 'Color', 'g', 'LineWidth', 2)
hold on;
line(logspace(0, log10(max(xvals)), 100), .1*logspace(0, log10(max(xvals)), 100), 'Color', 'r', 'LineWidth', 2)
errorbar(Ns, falseVals(:,1),  falseVals(:,1)-falseVals(:,2), falseVals(:,3)-falseVals(:,1), 'o');
set(gca, 'XScale', 'log')
ylim([0 ceil(max(falseVals(:,3))+1) ])
legend('One Percent', 'Ten Percent', 'Estimate with 95% CI', 'Location', 'NorthWest');
title(titleToUse, 'Interpreter', 'None');
xlabel('Input Mutation Count');
ylabel(ylabelToUse);
print('-dpng', [pngFilename '.png']);

close(gcf)