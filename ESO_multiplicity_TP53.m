ABS_maf=load_struct('/Users/amaro/Downloads/OLD_Eso_oncoprint/agg_maf');
absolute_table=load_struct('/Users/amaro/Downloads/OLD_Eso_Oncoprint/ESO_oldData_3.31.2014.ATW.ABSOLUTE.table.txt');
genes=load_struct('/Users/amaro/Downloads/cancer_gene_census.tsv');
all_genes=genes.gene;
absolute_table.Genomedoublings=str2double(absolute_table.Genomedoublings);

 for i=1:slength(ABS_maf)
                            strs=regexp(ABS_maf.sample{i},'-','split');
                            ABS_maf.sample{i}=strcat(strs{1},strs{2},'TPN');
 end
                        
                        
  for i=1:slength(absolute_table)
                         strs=regexp(absolute_table.array{i},'-','split');
                         absolute_table.array{i}=strcat(strs{1},strs{2},'TPN');
  end
  ABS_maf.modal_q_s=str2double(ABS_maf.modal_q_s);  
  ABS_maf=reorder_struct(ABS_maf,~ismember(ABS_maf.clonal_scna_mut_ix,'TRUE'));
  GDsamples=absolute_table.array(absolute_table.Genomedoublings>0);
  GD_maf=reorder_struct(ABS_maf,ismember(ABS_maf.sample,GDsamples));
  
%% mean plot  
figure
hold on
for i=1:length(all_genes)
k=ismember(GD_maf.Hugo_Symbol,all_genes{i});
if sum(k)>3
plot(i,mean(GD_maf.modal_q_s(k)),'k.');
text(i+2,mean(GD_maf.modal_q_s(k)),all_genes{i});
end
end

%% dist plots
% cgc genes
for i=1:length(all_genes)
    k=ismember(GD_maf.Hugo_Symbol,all_genes{i});
    if sum(k)>3
        dist_genes.(all_genes{i})=GD_maf.modal_q_s(k);
    end
end

distribution_ordered_plot(dist_genes.ARID1A,dist_genes.CDH11,dist_genes.CDKN2A,dist_genes.NOTCH2,dist_genes.NUP214,dist_genes.PDE4DIP...
    ,dist_genes.TP53,dist_genes.ZNF521)
ylabel('Mutliplicity')
axis([0 10 0 6])
xticklabel_rotate([1:8],45,fieldnames(dist_genes));

% more genes
all_genes=unique(GD_maf.Hugo_Symbol);
for i=1:length(all_genes)
    k=ismember(GD_maf.Hugo_Symbol,all_genes{i});
    if sum(k)>10
        dist_genes.(all_genes{i})=GD_maf.modal_q_s(k);
    end
end
fields=fieldnames(dist_genes)

distribution_ordered_plot(dist_genes.(fields{1}),dist_genes.(fields{2}),dist_genes.(fields{3}),dist_genes.(fields{4}),...
    dist_genes.(fields{5}),dist_genes.(fields{6}),dist_genes.(fields{7}),dist_genes.(fields{8}),dist_genes.(fields{9}),...
    dist_genes.(fields{10}),dist_genes.(fields{11}),dist_genes.(fields{12}),dist_genes.(fields{13}),dist_genes.(fields{14}),...
    dist_genes.(fields{15}),dist_genes.(fields{16}),dist_genes.(fields{17}),dist_genes.(fields{18}),dist_genes.(fields{19}),...
     dist_genes.(fields{20}),dist_genes.(fields{21}),dist_genes.(fields{22}),dist_genes.(fields{23}),dist_genes.(fields{24}))

xticklabel_rotate([1:24],45,fieldnames(dist_genes));


counts_mult = [0.080000000000000   0.750000000000000   0.166666666666667;...
   0.571428571428571   0.285714285714286   0.142857142857143];
   
   
P=bar(counts_mult,'stacked')
set(gca,'Ytick',[0,.5,1])
ylabel('Multiplicity of Mutations','FontSize',16)
set(gca,'XtickLabel',{'TP53','CDKN2A'});

set(P(1),'facecolor',[254/255 224/255 210/255],'EdgeColor',[1,1,1])
set(P(2),'facecolor',[252/255 146/255 114/255],'EdgeColor',[1,1,1])
set(P(3),'facecolor',[222/255 45/255 38/255],'EdgeColor',[1,1,1]) 
Pbaseline=get(P,'BaseLine');
set(Pbaseline{1}(1),'Color',[0 0 0],'LineWidth',.05);
box off
h=legend('Mult. 1','Mult. 2','Mult. 3plus','Location','NorthEastOutside');


