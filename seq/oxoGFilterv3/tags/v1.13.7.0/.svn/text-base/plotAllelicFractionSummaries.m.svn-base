function [] = plotAllelicFractionSummaries(outputFigureDir, titleString, basePngFilename, mafTable)
% [] = plotAllelicFractionSummaries(outputFigureDir, titleString, basePngFilename, mafTable)
%
%   outputFigureDir -- location to write png files.
%   titleString -- title of the generated plot
%   basePngFilename -- name of output png file.  Does not include the extension.  
%   mafTable -- tabular data as loaded by load_table.  This method looks
%       for the oxoGCut, i_t_ALT_F1R2, i_t_ALT_F1R2, i_t_REF_F1R2, i_t_REF_F1R2
%       columns.
% 
% See also load_table, myhistT, loadMAFTable

figure;
indices = logical(mafTable.oxoGCut);
alts = (mafTable.i_t_ALT_F1R2 + mafTable.i_t_ALT_F2R1);
af = alts ./ (alts + mafTable.i_t_REF_F1R2 + mafTable.i_t_REF_F2R1);

myhistT(af, .005:.01:.995, indices);
title(titleString, 'Interpreter', 'None')
xlabel('Allelic Fraction');
ylabel('Mutation Count');
legend('Filtered Mutations', 'Unfiltered Mutations', 'Location', 'North')
saveas(gcf, [outputFigureDir '/' basePngFilename '.png']);
close(gcf);