%% S. Figure 6a

seg1 = dataset('file','~/amaro/deTiN_Figure_Code_and_PDFs/Supplement/SupplementalFigure6/recapseg.BDD-CWGGV-Normal-SM-17M53.seg','ReadVarNames', true, 'Delimiter', '\t');
seg2 = dataset('file','~/amaro/deTiN_Figure_Code_and_PDFs/Supplement/SupplementalFigure6/recapseg.BDD-CWGGV-Tumor-SM-17M52.seg','ReadVarNames', true, 'Delimiter', '\t');

seg1.Segment_Mean = 2*2.^seg1.Segment_Mean;
seg2.Segment_Mean = 2*2.^seg2.Segment_Mean;

seg1.xstart = xhg19(seg_normal_brca.Chromosome,seg_normal_brca.Start);
seg1.xend = xhg19(seg_normal_brca.Chromosome,seg_normal_brca.End);
seg2.xstart = xhg19(seg_tumor_brca.Chromosome,seg_tumor_brca.Start);
seg2.xend = xhg19(seg_tumor_brca.Chromosome,seg_tumor_brca.End);


dx=.5;
figure()
hold on
load('rb_colomap.mat')

for j=1:length(seg1)
    
    x1=[ seg1.xstart(j) seg1.xstart(j) seg1.xend(j) seg1.xend(j) seg1.xstart(j)];
    y1=[1-dx 1+dx 1+dx 1-dx 1-dx];
    
    if seg1.Segment_Mean(j) > 4
        
        c=[256,0,0];
    else
        c=rb_colormap(round((seg1.Segment_Mean(j)/4)*255)+1,:);
    end
   
    
    patch(x1,y1,c./256,'edgecolor','none')
end
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[.5,1.5],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
set(gca,'ytick',[])

title('BDD-CWGGV-Normal')

figure()
hold on
for j=1:length(seg2)
    
    x1=[ seg2.xstart(j) seg2.xstart(j) seg2.xend(j) seg2.xend(j) seg2.xstart(j)];
    y1=[1-dx 1+dx 1+dx 1-dx 1-dx];
    
    if seg2.Segment_Mean(j) > 4
        
        c=[256,0,0];
    else
        c=rb_colormap(round((seg2.Segment_Mean(j)/4)*255)+1,:);
    end
   
    
    patch(x1,y1,c./256,'edgecolor','none')
end
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[.5,1.5],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
set(gca,'ytick',[])
title('BDD-CWGGV-Tumor')

start_snp=126853103;	end_snp=262093187;
SNPs = readtable('~/amaro/deTiN_Figure_Code_and_PDFs/Supplement/SupplementalFigure6/BDD-CWGGV-TP-NT-SM-17M52-SM-17M53.call_stats.txt','Delimiter','\t','HeaderLines',1);
SNPs.contig = chromosome2num_legacy(SNPs.contig);
SNPs = SNPs(SNPs.position> start_snp & SNPs.position<end_snp & SNPs.contig ==1 & SNPs.tumor_f > .2& SNPs.tumor_f < .9,:);
SNPs.direction = SNPs.normal_f>.5;
hold on
plot(SNPs.position(SNPs.direction==1),SNPs.normal_f(SNPs.direction==1),'r.','MarkerSize',10)
plot(SNPs.position(SNPs.direction==0),SNPs.normal_f(SNPs.direction==0),'b.','MarkerSize',10)
figure()
hold on
plot(SNPs.position(SNPs.direction==1),SNPs.tumor_f(SNPs.direction==1),'r.','MarkerSize',10)
plot(SNPs.position(SNPs.direction==0),SNPs.tumor_f(SNPs.direction==0),'b.','MarkerSize',10)
%% S. Figure 6b