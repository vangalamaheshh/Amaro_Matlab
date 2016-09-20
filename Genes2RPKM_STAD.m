sig_genes=load_struct('/xchip/cga/gdac-prod/cga/jobResults/MutSig_S2N/Normal_v4/3851063/Normal_v4.sig_genes.txt')
% RPKM_mat=load_struct('~/Projects/STAD/STAD_RPKM_vals.txt')
%  fields=fieldnames(RPKM_mat);
% 
%  for i=1:slength(RPKM_mat)
%  for j=2:length(fields)-3
%  vals(j)=str2num(RPKM_mat.(fields{j}){i});
%  end
%  RPKM_stats.median{i,1}=nanmedian(vals);  
%  RPKM_stats.nanstd{i,1}=std(vals);        
%  RPKM_stats.gene{i,1}=RPKM_mat.gene{i};     
%  end
% save('~/Projects/STAD/RPKM_stats.mat',RPKM_stats');

load('~/Projects/STAD/RPKM_stats.mat')
sig_genes.q=cellfun(@str2num,sig_genes.q); 
sig_genes=reorder_struct(sig_genes,(sig_genes.q<.05)); 

RPKM_stats.median=cell2mat(RPKM_stats.median); 
for i=1:slength(sig_genes)
loc=find(ismember(RPKM_stats.gene,sig_genes.gene{i})==1);
sig_genes.median_RPKM(i,1)=mean(RPKM_stats.median(loc));
end

