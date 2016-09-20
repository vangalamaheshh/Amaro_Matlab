function GCT_check_same_allele(seg1,seg2,cov1,cov2)

seg1=rmfield(seg1,{'headline','header'});
seg2=rmfield(seg2,{'headline','header'});
if ~isnumeric(seg1.Chromosome)
seg1.Chromosome=chromosome2num_legacy(seg1.Chromosome);
end
if ~isnumeric(seg2.Chromosome)
seg2.Chromosome=chromosome2num_legacy(seg2.Chromosome);
end
seg1.xstart=xhg19(seg1.Chromosome,seg1.Start_bp);
seg2.xstart=xhg19(seg2.Chromosome,seg2.Start_bp);
seg1.xend=xhg19(seg1.Chromosome,seg1.End_bp);
seg2.xend=xhg19(seg2.Chromosome,seg2.End_bp);

cov1=rmfields_if_exist(cov1,{'headline','header'});
cov2=rmfields_if_exist(cov2,{'headline','header'});
cov1.Chromosome=chromosome2num_legacy(cov1.Chromosome);
cov2.Chromosome=chromosome2num_legacy(cov2.Chromosome);
seg1=reorder_struct(seg1,seg1.n_hets>20);

for seg=1:slength(seg1)
c1=reorder_struct(cov1,cov1.Start_position>seg1.Start_bp(seg)&cov1.Start_position<seg1.End_bp(seg)&cov1.Chromosome==seg1.Chromosome(seg));
c2=reorder_struct(cov2,cov2.Start_position>seg1.Start_bp(seg)&cov2.Start_position<seg1.End_bp(seg)&cov2.Chromosome==seg1.Chromosome(seg));
c1=reorder_struct(c1,ismember(c1.Start_position,c2.Start_position));
c2=reorder_struct(c2,ismember(c2.Start_position,c1.Start_position));
cshared(seg,1)=min(min(corrcoef(c1.i_t_alt_count./(c1.i_t_alt_count+c1.i_t_ref_count),c2.i_t_alt_count./(c2.i_t_alt_count+c2.i_t_ref_count))));
end



load('rb_colomap.mat')

figure()
hold on
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[-1,1],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end

for i=1:slength(seg1)
plot([seg1.xstart(i),seg1.xend(i)],[cshared(i),cshared(i)],'b-')
end
end

function example
seg1=load_table('~/Downloads/TCG-Testes_DFCI_9-Tumor-SM-4PDDM.tsv');
seg2=load_table('~/Downloads/TCG-Testes_DFCI_9-Tumor-SM-4PDDN.tsv');
cov1=load_table('~/Downloads/4PDDM.cov');
cov2=load_table('~/Downloads/4PDDN.cov');
GCT_check_same_allele(seg1,seg2,cov1,cov2)

seg1=load_table('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDF.tsv');
seg2=load_table('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDH.tsv');
cov1=load_table('~/Downloads/4PDDF.cov');
cov2=load_table('~/Downloads/4PDDH.cov');
GCT_check_same_allele(seg1,seg2,cov1,cov2)


seg1=load_table('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDF.tsv');
seg2=load_table('~/Downloads/TCG-Testes_DFCI_7-Tumor-SM-4PDDG.tsv');
cov1=load_table('~/Downloads/4PDDF.cov');
cov2=load_table('~/Downloads/4PDDG.cov');
GCT_check_same_allele(seg1,seg2,cov1,cov2)
end