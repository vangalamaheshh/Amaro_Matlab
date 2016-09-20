x=load_struct('/xchip/cga_home/amaro/Barretts/ABSOLUTE_Run_Paper_Combined_Sets/MafFiles_for_figures_Barretts_Pairs.txt');

M=load_struct(x.maf{1});

for i=2:slength(x)
m=load_struct(x.maf{i});
M=mergeStruct(M,m);
M=rmfield(M,'N');
if mod(i,5)==0
    i
end
end
    

save_struct(M,'/Users/amaro/Documents/BarrettsRevisionFigures/ABSOLUTE.merged.maf')
sig_genes=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/genes_for_ESO.txt');
OncoM=reorder_struct(M,ismember(M.Hugo_Symbol,sig_genes.gene));
save_struct(M,'/Users/amaro/Documents/BarrettsRevisionFigures/ABSOLUTE.OncoGenes.maf')



clear M m x

x=load_struct('/xchip/cga_home/amaro/CLL/ABSOLUTE/TimePoint2.mafs');

M=load_struct(x.maf{1});

for i=2:slength(x)
m=load_struct(x.maf{i});
M=mergeStruct(M,m);
M=rmfield(M,'N');
if mod(i,5)==0
    i
end
end
    

save_struct(M,'/xchip/cga_home/amaro/CLL/ABSOLUTE/CombinedABSResultsForDanTPLater.maf')
