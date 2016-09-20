M=load_struct('Filtered_GCT_UK.maf');
samples=unique(M.Tumor_Sample_Barcode);
hist_table=load_struct('~/Documents/UK_histology_table.txt');
ABS_SEG=load_struct('/Volumes/xchip_cga_home/amaro/ForEli/FreezeSet/UK_ABS_SEGs/AggregatedUKABS.seg'); ABS_SEG=rmfields_if_exist(ABS_SEG,{'header','headline'});
ABS_SEG.Chromosome=chromosome2num_legacy(ABS_SEG.Chromosome);
ABS_SEG.Startbp=str2double(ABS_SEG.Startbp);ABS_SEG.Endbp=str2double(ABS_SEG.Endbp);
Events = {'KRAS','RPL5','KIT','APC','12p','KRAS_amp'};
ABS_SEG=reorder_struct(ABS_SEG,~ismember(ABS_SEG.sample,'UK_124-Tumor2-TP-NB'));
for i=1:slength(M)
    M.sample{i,1}=regexp(M.Tumor_Sample_Barcode{i},'[0-9]+','match');
    M.sample{i,1}=char(M.sample{i});
end

for i=1:slength(ABS_SEG)
    ABS_SEG.sample_num{i,1}=regexp(ABS_SEG.sample{i},'[0-9]+','match');
    ABS_SEG.sample_num{i,1}=char(ABS_SEG.sample_num{i});
end
CHR12_SEG=reorder_struct(ABS_SEG,ABS_SEG.Chromosome==12&ABS_SEG.Startbp<34852809);
KRAS_SEG=reorder_struct(ABS_SEG,ABS_SEG.Chromosome==12&ABS_SEG.Startbp<25356180&ABS_SEG.Endbp>25356180);

M=reorder_struct(M,~ismember(M.Tumor_Sample_Barcode,{'UK_124-Tumor2','UK_111-Tumor2'}));
hist_table=reorder_struct(hist_table,ismember(hist_table.Tumor_Number,ABS_SEG.sample_num));
M=reorder_struct(M,ismember(M.sample,ABS_SEG.sample_num));

Events = {'KRAS','RPL5','KIT','APC','12p','KRAS_amp'};
load('rb_colomap.mat')
colors=[236,250,80;80,236,250;250,5,128;5,250,127];



dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;

N=length(hist_table.Tumor_Number);
NV=length(Events);



missense=[179/255 222/255 105/255];
frameshift=[102/255 0/255 204/255];
nonsense=[153/255 0 13/255];
splice=[1 1 51/255];
silent=[237/255 248/255 251/255];
hold on
set(gca,'visible','off');


for v=1:NV
    for s=1:N
        %these correspond to the shapes of the panels in this case rectangles
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
        
        patch(x1,y1,[.95 .95 .95],'edgecolor',[.95 .95 .95])
        
    end
end

for v=1:NV
    event=Events{v};
    for s=1:N
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
        sample=hist_table.Tumor_Number{s};
        if ismember(sample,M.sample(ismember(M.Hugo_Symbol,event)))
            type=M.Variant_Classification(ismember(M.sample,sample)&ismember(M.Hugo_Symbol,event));
            if ismember('Nonsense_Mutation',type)
                patch(x1,y1,nonsense,'edgecolor',nonsense)
            elseif ismember('Splice_Site',type)
                patch(x1,y1,splice,'edgecolor',splice)
                
            elseif ismember('Frame_Shift_Del',type)||ismember('Frame_Shift_Ins',type)
                patch(x1,y1,frameshift,'edgecolor',frameshift)
            elseif ismember('Missense_Mutation',type)||ismember('Nonstop_Mutation',type)
                patch(x1,y1,missense,'edgecolor',missense)
            elseif ismember('Silent',type)||ismember('3''UTR',type)||ismember('5''UTR',type)
                patch(x1,y1,silent,'edgecolor',silent)
            elseif ismember('In_Frame_Del',type)||ismember('In_Frame_Ins',type)
                patch(x1,y1,inframe,'edgecolor',inframe)
            end
        end
        if isequal(event,'12p')
            seg=reorder_struct(CHR12_SEG,ismember(CHR12_SEG.sample_num,sample));
            seg.modala1=str2double(seg.modala1); seg.modala2=str2double(seg.modala2);
            tau=nanmedian(seg.modala1+seg.modala2); 
            if tau>8 
                tau=8;
            end
            if tau >= 2.5
                c=rb_colormap(round((tau/8)*127)+127,:);
            elseif tau <= 2.5
                
            c=[255 255 255];
            end
            patch(x1,y1,c./255,'edgecolor','none')
        end
        if isequal(event,'KRAS_amp')
            seg=reorder_struct(CHR12_SEG,ismember(CHR12_SEG.sample_num,sample)); seg.rescaled_total_cn=str2double(seg.rescaled_total_cn);
            tau=nanmedian(seg.rescaled_total_cn);
            kras_cn=reorder_struct(KRAS_SEG,ismember(KRAS_SEG.sample_num,sample)); kras_cn.rescaled_total_cn=str2double(kras_cn.rescaled_total_cn);
            if slength(kras_cn)>0
                
                cn=max(kras_cn.rescaled_total_cn);
               % [cn,tau]
                if cn/tau>=2 || cn > 10
                    patch(x1,y1,[1 0 0],'edgecolor','none')
                end
            else
                sample
                patch(x1,y1,[.5 .5 .5],'edgecolor','none')
            end
        end
        
    end
end

extra_tracks=1;
hist_table.HistClass(ismember(hist_table.Histology,'Seminoma'))={1};
hist_table.HistClass(ismember(hist_table.Histology,'Non-Seminoma'))={2};
hist_table.HistClass(ismember(hist_table.Histology,'Mixed'))={3};
hist_table.HistClass(ismember(hist_table.Histology,'None'))={4};
for v=1:extra_tracks;
    for s=1:N
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v+extra_tracks;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    if v==1
        for z=1:3
        if ismember(z,hist_table.HistClass{s})
            patch(x1,y1,colors(z,:)./255,'edgecolor','none')
        else
            if sum(ismember(hist_table.HistClass{s},[1,2,3]))==0
                patch(x1,y1,colors(4,:)./255,'edgecolor','none')
            end
        end
        end
    end
    end
end
    
    