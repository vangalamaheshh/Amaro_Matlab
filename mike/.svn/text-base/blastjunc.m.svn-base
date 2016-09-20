function [B,F]=blastjunc(id,parameters)
if ~exist('parameters','var'), parameters='-e1e-50 -FF'; end
F=load_fasta(['/xchip/tcga/gbm/analysis/lawrence/wgs/jump/gg/cp/' id '/ClosePairs.assembly.brief.fasta']);
B=blast(F.seq{1},'hg18',parameters);
