SEG=load_struct('/Volumes/xchip_cga_home/amaro/ForEli/GCT_ACS/PrimaryMutSetGCT.seg')
for i=1:slength(SEG)
    SEG.sample{i,1}=regexp(SEG.sample{i},'DFCI_[0-9]+','match');
    SEG.sample{i,1}=char(SEG.sample{i});
end
ClinicalData=load_struct('GCT_Armlvl_data/ClinicalPathData.txt');
for i=1:slength(ClinicalData)
ClinicalData.sample{i,1}=regexp(ClinicalData.FH_names{i},'DFCI_[0-9]+','match');
ClinicalData.sample{i,1}=char(ClinicalData.sample{i});
end

ClinicalData=reorder_struct(ClinicalData,ismember(ClinicalData.sample,SEG.sample));


for i=1:slength(ClinicalData)
    
    ClinicalData.HistClass{i,:}=str2double(regexp(ClinicalData.OrchHistology{i},'[0-9]+','match'));
    
end



SEG=rmfield_if_exist(SEG,'headline');
SEG=rmfield_if_exist(SEG,'header');
SEG.Chromosome=chromosome2num_legacy(SEG.Chromosome);
SEG.Startbp=str2double(SEG.Startbp);
SEG.Endbp=str2double(SEG.Endbp);
SEG.xStart=xhg19(SEG.Chromosome,SEG.Startbp);
SEG.tau=str2double(SEG.tau);
SEG.f=str2double(SEG.f);
SEG.t_n=SEG.tau;
SEG.t_n(SEG.tau>4)=4;
SEG.n_hets=str2double(SEG.n_hets);
SEG.n_probes=str2double(SEG.n_probes);
%SEG=reorder_struct(SEG,SEG.n_hets>1);
SEG=reorder_struct(SEG,~isnan(SEG.f));
SEG=reorder_struct(SEG,SEG.n_probes>10);
SEG=reorder_struct(SEG,SEG.Chromosome<23);
load('rb_colomap.mat')
load('grncolormap.mat')
SEG.xEnd=xhg19(SEG.Chromosome,SEG.Endbp);
starts=xhg19(1:22,0*(1:22),'hg19');
Is=ClinicalData.sample;
figure()
%VHL=xhg19([3 3],[10182692 10193904]);
hold on
% for i=1:length(starts)
% 
%     plot([1 length(Is)],[starts(i) starts(i)],'linestyle','-','color',0.5*[1 1 1])
%     
% end
%    


% for i=1:length(Is)+1
%     plot([i i],[0 max(SEG.xEnd)],'linestyle','--','color',0.5*[1 1 1])
% end
dx=.5;
for i=2:2:22
    s=.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    if i<22
    y1=[starts(i) starts(i) starts(i+1) starts(i+1) starts(i)];
    else
        y1=[starts(i) starts(i) SEG.xEnd(end) SEG.xEnd(end) starts(i)];
    end
    patch(x1,y1,'k','edgecolor','none')   
end


for i=1:2:22
    s=.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    if i<22
    y1=[starts(i) starts(i) starts(i+1) starts(i+1) starts(i)];
    else
        y1=[starts(i) starts(i) SEG.xEnd(end) SEG.xEnd(end) starts(i)];
    end
    patch(x1,y1,[.5 .5 .5],'edgecolor','none')   
end

for i=1:length(Is)  
    I_seg=reorder_struct(SEG,ismember(SEG.sample,Is{i}));
    s=i+.5;
    for j=1:slength(I_seg)
        y1=[ I_seg.xStart(j) I_seg.xStart(j) I_seg.xEnd(j) I_seg.xEnd(j) I_seg.xStart(j)];
        x1=[s-dx s+dx s+dx s-dx s-dx];
        if I_seg.t_n(j)>1.9 && I_seg.t_n(j)<2.1 && I_seg.f(j)<.46
                    c=grncolor(round((1-(I_seg.f(j)/.5))*255)+1,:);
           patch(x1,y1,c./255,'edgecolor','none') 
        else

           c=rb_colormap(round((I_seg.t_n(j)/4)*255)+1,:);
         patch(x1,y1,c./255,'edgecolor','none')    
        end
          
       
    end
end


for i=1:length(starts)
    if i<22
l(i)=((starts(i+1)-starts(i))/2)+starts(i);
    else
        l(i)=starts(i)+((SEG.xEnd(end)-starts(i))/2);
    end
end

set(gca,'ytick',l(2:2:22),'ytickLabel',num2str([2:2:22]'));


extra_track_start=max(I_seg.xEnd)+100000;
colors=[236,250,80;80,236,250;250,5,128;5,250,127];
for i=1:slength(ClinicalData)
    s=i+.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    y1=[ extra_track_start extra_track_start extra_track_start+100000000 extra_track_start+100000000 extra_track_start];
    patch(x1,y1,[.5 .5 .5],'edgecolor','none')
    for j=1:3
        if ismember(j,ClinicalData.HistClass{i})
            patch(x1,y1,colors(j,:)./255,'edgecolor','none')
        else
            if sum(ismember(ClinicalData.HistClass{i},[1,2,3]))==0
                patch(x1,y1,colors(4,:)./255,'edgecolor','none')
            end
        end
    end
end

extra_track_start=extra_track_start+100000000;
for i=1:slength(ClinicalData)
    s=i+.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    y1=[ extra_track_start extra_track_start extra_track_start+100000000 extra_track_start+100000000 extra_track_start];
    patch(x1,y1,[.5 .5 .5],'edgecolor','none')
    if isequal(ClinicalData.Vitalstatus{i},'1')
        patch(x1,y1,[1 0 0],'edgecolor','none')
    else
        patch(x1,y1,[0 0 0],'edgecolor','none')
    end
end




%set(gca,'xtick',[1.5:1:35.5])
box off
xlim([0 length(Is)+1])
ylim([-1 extra_track_start+100000000])
%  plot([1 11],[VHL(1) VHL(1)],'linestyle','-','color',0*[1 1 0])
%     plot([1 11],[VHL(2) VHL(2)],'linestyle','-','color',0*[1 1 0])

xticklabel_rotate_simplified([1:length(unique(SEG.sample))]+.5,90,Is)
