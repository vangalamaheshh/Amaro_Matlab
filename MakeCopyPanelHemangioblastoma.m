SEG=load_table('/Users/amaro/Downloads/AggregateAllelicSegments.seg');
SEG=rmfield_if_exist(SEG,'headline');
SEG=rmfield_if_exist(SEG,'header');
SEG.xStart=xhg19(SEG.Chromosome,SEG.Start_bp);
SEG.xEnd=xhg19(SEG.Chromosome,SEG.End_bp);
starts=xhg19(1:24,0*(1:24),'hg19');
Is=unique(SEG.Sample);
figure()
VHL=xhg19([3 3],[10182692 10193904]);
hold on
for i=1:length(starts)

    plot([1 11],[starts(i) starts(i)],'linestyle','-','color',0.5*[1 1 1])
    
end
   


for i=1:length(Is)+1
    plot([i i],[0 max(SEG.xEnd)],'linestyle','--','color',0.5*[1 1 1])
end
dx=.5;
for i=1:2:24
    s=.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    y1=[starts(i) starts(i) starts(i+1) starts(i+1) starts(i)];
    patch(x1,y1,'k','edgecolor','none')   
end

for i=1:length(Is)  
    I_seg=reorder_struct(SEG,ismember(SEG.Sample,Is{i}));
    s=i+.5;
    for j=1:slength(I_seg)
    if I_seg.allelic_copy_ratio(j)<1.85 
        y1=[ I_seg.xStart(j) I_seg.xStart(j) I_seg.xEnd(j) I_seg.xEnd(j) I_seg.xStart(j)];
        x1=[s-dx s+dx s+dx s-dx s-dx];
        if I_seg.allelic_copy_ratio(j)<1.5
            patch(x1,y1,[0 0 1],'edgecolor','none')
            
        else
        c=255*(.5-(2-I_seg.allelic_copy_ratio(j)))/255;
        patch(x1,y1,[c c 255/255],'edgecolor','none')   
        end
    elseif I_seg.allelic_copy_ratio(j)>2.15
        y1=[ I_seg.xStart(j) I_seg.xStart(j) I_seg.xEnd(j) I_seg.xEnd(j) I_seg.xStart(j)];
        x1=[s-dx s+dx s+dx s-dx s-dx];
        if I_seg.allelic_copy_ratio(j)>2.5
                        patch(x1,y1,[1 0 0],'edgecolor','none')
        else
            c=255*(.5-(I_seg.allelic_copy_ratio(j)-2))/255;
        patch(x1,y1,[255/255 c c],'edgecolor','none')   
        end

    
    end
    end
end
for i=1:length(starts)-1
l(i)=((starts(i+1)-starts(i))/2)+starts(i);
end
set(gca,'ytick',l(1:2:23),'ytickLabel',num2str([1:2:23]'));
set(gca,'xtick',[])
box off
xlim([0 11])
ylim([-1 max(starts)])
 plot([1 11],[VHL(1) VHL(1)],'linestyle','-','color',0*[1 1 0])
    plot([1 11],[VHL(2) VHL(2)],'linestyle','-','color',0*[1 1 0])


