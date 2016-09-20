function MakeCopyPanelFromAllelicseg(SEG)
SEG=rmfield_if_exist(SEG,'headline');
SEG=rmfield_if_exist(SEG,'header');
SEG.xStart=xhg19(SEG.Chromosome,SEG.Start);

SEG.Chromosome=chromosome2num_legacy(SEG.Chromosome);
SEG.Segment_Mean=str2double(SEG.Segment_Mean);
SEG.t_n=SEG.Segment_Mean;
SEG.t_n(SEG.Segment_Mean>1.5)=1.5;
SEG.t_n(SEG.Segment_Mean<-1.5)=-1.5;
SEG.t_n=SEG.t_n+1.5;
SEG=reorder_struct(SEG,SEG.Chromosome<23);
load('rb_colomap.mat')
SEG.xEnd=xhg19(SEG.Chromosome,SEG.End);
starts=xhg19(1:22,0*(1:22),'hg19');
Is=unique(SEG.Sample,'stable');
figure()
%VHL=xhg19([3 3],[10182692 10193904]);
hold on
% for i=1:length(starts)
% 
%     plot([1 length(Is)],[starts(i) starts(i)],'linestyle','-','color',0.5*[1 1 1])
%     
% end
%    

% 
% for i=1:length(Is)+1
%     plot([i i],[0 max(SEG.xEnd)],'linestyle','--','color',0.5*[1 1 1])
% end
samples=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/SamplesTableWithClusteringOrder.txt');
Is=samples.pair_id;
dx=.5;
%%%%%%%%%%%% Code for Barretts Project in particular
barretts=[218/255 112/255 214/255];
eso=[34/255 139/255 34/255];
ND=[255/255 182/255 193/255];

samples.pair_id=tr(samples.pair_id,'.','-');
for i=1:slength(samples)
    %k=find(ismember(samples.pair_id,Is{i}));
    sampletype=samples.tissue{i};
    doubled=samples.GD{i};
    s=i+.5;
    y1=[-10^8*2 -10^8*2 -10^8 -10^8 -10^8*2];
    x1=[s-dx s+dx s+dx s-dx s-dx];
    I_seg=reorder_struct(SEG,ismember(SEG.Sample,samples.pair_id{i}));
    if isequal(sampletype,'ESO')
        patch(x1,y1,eso,'edgecolor','none')
    elseif isequal(sampletype,'ND')
        patch(x1,y1,ND,'edgecolor','none')
    else
       patch(x1,y1,barretts,'edgecolor','none')

    end
    y1=[-10^8 -10^8 -1 -1 -10^8];
    if isequal(doubled,'0')
    else
      patch(x1,y1,[0 0 0],'edgecolor','none')
    end
         
end



for i=2:2:22
    s=.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    if i<22
    y1=[starts(i) starts(i) starts(i+1) starts(i+1) starts(i)];
    else
        y1=[starts(i) starts(i) max(SEG.xEnd) max(SEG.xEnd) starts(i)];
    end
    patch(x1,y1,'k','edgecolor','none')   
end
for i=1:2:22
    s=.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    if i<22
    y1=[starts(i) starts(i) starts(i+1) starts(i+1) starts(i)];
    else
        y1=[starts(i) starts(i) max(SEG.xEnd) max(SEG.xEnd) starts(i)];
    end
    patch(x1,y1,[.5 .5 .5],'edgecolor','none')   
end


for i=1:slength(samples)  
    I_seg=reorder_struct(SEG,ismember(SEG.Sample,samples.pair_id{i}));
    s=i+.5;
    for j=1:slength(I_seg)
        y1=[ I_seg.xStart(j) I_seg.xStart(j) I_seg.xEnd(j) I_seg.xEnd(j) I_seg.xStart(j)];
        x1=[s-dx s+dx s+dx s-dx s-dx];
        c=rb_colormap(round(((I_seg.t_n(j)/3))*255)+1,:);
        if I_seg.Segment_Mean(j)<.1&&I_seg.Segment_Mean(j)>-.1
            c=[255 255 255];
        end
        
        patch(x1,y1,c./255,'edgecolor','none')   
       
    end
end


for i=1:length(starts)
    if i<22
l(i)=starts(i)+((starts(i+1)-starts(i))/2);
    else
        l(i)=starts(i)+(max(SEG.xEnd)-starts(i))/2;
    end
end



%  plot([1 11],[VHL(1) VHL(1)],'linestyle','-','color',0*[1 1 0])
%     plot([1 11],[VHL(2) VHL(2)],'linestyle','-','color',0*[1 1 0])



strs=char2cell(num2str([1:22]'));
set(gca,'ytick',[-10^8*1.5,-10^8*.5,l(1:22)],'ytickLabel',['Sample Type';'GD';strs]);
set(gca,'xtick',[1.5:1:length(Is)])
%set(gca,'xtick',[])
for i=1:length(Is)
strs=split(Is{i},'.');
ind{i,1}=strs{1};
end
box off
xlim([0 length(Is)+1])
ylim([-10^8*2 max(SEG.xEnd)])



end


function test

barretts=[218/255 112/255 214/255];
eso=[34/255 139/255 34/255];



end

