CNV_HESC=load_table('~/Downloads/CNVs_suppTable3_4-2.txt');
CNV_HESC=rmfield(CNV_HESC,{'headline','header'});
samples=unique(CNV_HESC.Sample_ID);
CNV_HESC.xstart=xhg19(chromosome2num_legacy(CNV_HESC.Chromosome),CNV_HESC.Start_bp);
CNV_HESC.xend=xhg19(chromosome2num_legacy(CNV_HESC.Chromosome),CNV_HESC.End_bp);
load('rb_colomap.mat')

dx=.5;
%figure()

figure()
hold on
for i=1:length(samples)
    seg1=reorder_struct(CNV_HESC,ismember(CNV_HESC.Sample_ID,samples{i}));
    for j=1:slength(seg1)
        x1=[ seg1.xstart(j) seg1.xstart(j) seg1.xend(j) seg1.xend(j) seg1.xstart(j)];
        y1=[i-dx i+dx i+dx i-dx i-dx];
        c=rb_colormap(round(seg1.Copy_Number(j)/4)*255+1,:);
        if seg1.Copy_Number(j)<2
            c=[0 0 256];
        else
            c=[256 0 0];
        end
        patch(x1,y1,c./256,'edgecolor','none')
        
    end
end


aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[0,38],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
set(gca,'ytick',[])