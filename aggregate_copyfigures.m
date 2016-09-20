load('~/Documents/final_code_folder_GCT/FigureWorkSpaceForCN.mat')

for i=1:22
violin_cell{i,1}=violins_amp_full(violins_count==i)+1;
end

violin({violin_cell{2:end}},'mc',[],'medc','k')
set(gca,'XTickLabel',{xlabs1{2:end}});
set(gca,'XTick',1:length(xlabs1)-1);
ylabel({'# of arms with at least 1 allele amplified','(Excluding WGD events)'})
