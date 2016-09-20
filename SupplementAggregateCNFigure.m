ABS_SEG=load_table('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/ABS_Segs/AggregatedGCTABS.seg'); ABS_SEG=rmfield(ABS_SEG,{'header','headline'});
ABS_Table=load_table('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/AggregatedGCTABSTable.txt')
load('rb_colomap.mat')

samples=unique(ABS_SEG.sample);
index = 1;
figure()
hold on
for i = 1:length(samples)
    
    samp1=samples{i};
    seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
    GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
    absolute_allelic_cn_plot( seg1 , GD , index)
    index = index + 2;
    
end
index = 1;
for i = 1:length(samples)
    line([1,2881033286.00000],[index+.5,index+.5],'color',[0 0 0],'linestyle','--','LineWidth',.5)
    index = index + 2;
end
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[0,121],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end


