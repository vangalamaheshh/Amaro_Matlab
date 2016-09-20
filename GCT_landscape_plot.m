M=load_struct('~/GitHub/Testes_Submission/GCT_Final_Tables_Figures/Supplementary_Table_3.maf');
M_ABS=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/ABS_Mafs/AggregatedGCT.maf');

SEG=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/GCT_ACS/PrimaryMutSetGCT.seg');
ABS_SEG=load_table('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/ABS_Segs/AggregatedGCTABS.seg'); ABS_SEG=rmfield(ABS_SEG,{'header','headline'});
KRAS.Chr=12; KRAS.start=25356180; KRAS.end=25405854;
MSIG_Pairs=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/MsigPairtable.tsv');
M.Tumor_Sample_Barcode=M.sample;
M=reorder_struct(M,ismember(M.Tumor_Sample_Barcode,MSIG_Pairs.case_sample));
M_ABS=reorder_struct(M_ABS,ismember(M_ABS.sample_name,MSIG_Pairs.case_sample));
ABS_Table=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/AggregatedGCTABSTable.txt');
ClinicalData=load_struct('~/Documents/GCT_Final_STable1.txt');

for i=1:slength(SEG)
    SEG.sample{i,1}=regexp(SEG.sample{i},'DFCI_[0-9]+','match');
    SEG.sample{i,1}=char(SEG.sample{i}); 
end
ABS_Table=reorder_struct(ABS_Table,ismember(ABS_Table.array,MSIG_Pairs.pair_id));
SEG.Chromosome=chromosome2num_legacy(SEG.Chromosome);
SEG.Startbp=str2double(SEG.Startbp);
SEG.Endbp=str2double(SEG.Endbp);
SEG.xStart=xhg19(SEG.Chromosome,SEG.Startbp);
SEG.tau=str2double(SEG.tau);
SEG.t_n=SEG.tau;
SEG.t_n(SEG.tau>4)=4;

ClinicalData.sample=ClinicalData.patient_id;
ABS_SEG.pair_id=ABS_SEG.sample;
for i=1:slength(ABS_SEG)
    ABS_SEG.sample{i,1}=regexp(ABS_SEG.pair_id{i},'DFCI_[0-9]+','match');
    ABS_SEG.sample{i,1}=char(ABS_SEG.sample{i});
end
CHR12_SEG=reorder_struct(ABS_SEG,ABS_SEG.Chromosome==12&ABS_SEG.Startbp<34852809);
KRAS_SEG=reorder_struct(ABS_SEG,ABS_SEG.Chromosome==12&ABS_SEG.Startbp<25356180&ABS_SEG.Endbp>25356180);
M.patient=M.sample;
for i=1:slength(M)
    M.sample{i,1}=regexp(M.patient{i},'DFCI_[0-9]+','match');
    M.sample{i,1}=char(M.sample{i});
end
for i=1:slength(M_ABS)
    M_ABS.sample{i,1}=regexp(M_ABS.sample{i},'DFCI_[0-9]+','match');
    M_ABS.sample{i,1}=char(M_ABS.sample{i});
end

for i=1:slength(ABS_Table)
    ABS_Table.sample{i,1}=regexp(ABS_Table.array{i},'DFCI_[0-9]+','match');
    ABS_Table.sample{i,1}=char(ABS_Table.sample{i});
end

pset_list=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/JointMaf/pset_list.tsv');
M=reorder_struct(M,ismember(M.patient,pset_list.case_sample));
M_ABS=reorder_struct(M_ABS,ismember(M_ABS.sample,M.sample));
ClinicalData=reorder_struct(ClinicalData,ismember(ClinicalData.sample,M.sample));
SEG=reorder_struct(SEG,ismember(SEG.sample,M.sample));
KRAS_SEG=reorder_struct(KRAS_SEG,ismember(KRAS_SEG.sample,M.sample));
Events = {'KRAS','RPL5','TP53','12p','KRAS_amp'};
load('rb_colomap.mat')
colors=[236,250,80;80,236,250;250,5,128;5,250,127];
figure()
hold on
missense=[179/255 222/255 105/255];
frameshift=[102/255 0/255 204/255];
nonsense=[153/255 0 13/255];
splice=[1 1 51/255];
silent=[237/255 248/255 251/255];
inframe=[0 0 1];
hold on
set(gca,'visible','off');
%setting dims of boxs
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;

N=length(unique(ClinicalData.sample));
NV=length(Events);

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
        sample=ClinicalData.sample{s};
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
            seg=reorder_struct(CHR12_SEG,ismember(CHR12_SEG.sample,sample));
           
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
            seg=reorder_struct(CHR12_SEG,ismember(CHR12_SEG.sample,sample));
            tau=median(seg.rescaled_total_cn);
            kras_cn=reorder_struct(KRAS_SEG,ismember(KRAS_SEG.sample,sample));
            if slength(kras_cn)>0
                
                cn=max(kras_cn.rescaled_total_cn);
               % [cn,tau]
                if cn/tau>=2 
                    patch(x1,y1,[1 0 0],'edgecolor','none')
                end
            else
                sample
                patch(x1,y1,[.5 .5 .5],'edgecolor','none')
            end
        end
        
    end
end
extra_tracks=4;
ABS_Table.Genomedoublings=str2double(ABS_Table.Genomedoublings);
for v=1:extra_tracks;
    for s=1:N
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v+extra_tracks;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    patch(x1,y1,[.5 .5 .5],'edgecolor','none')
    if isequal(ClinicalData.Vitalstatus{s},'2') && v==2
        patch(x1,y1,[1 0 0],'edgecolor','none')
    elseif v==1 && isequal(ClinicalData.Vitalstatus{s},'1')
        patch(x1,y1,[0 0 0],'edgecolor','none')
    end
    if v==1
        if isequal(ClinicalData.Resistance{s},'0')
            patch(x1,y1,[0.5 0.5 0.5],'edgecolor','none')
        elseif isequal(ClinicalData.Resistance{s},'1')
            patch(x1,y1,[0 0 0],'edgecolor','none')
        end
    end
        
    if v==3
            if isequal(ClinicalData.Histology{s},'Seminoma')
            patch(x1,y1,[255,255,102]./255,'edgecolor','none')
            elseif isequal(ClinicalData.Histology{s},'Non-seminoma - mixed')
                patch(x1,y1,[102,204,204]./255,'edgecolor','none')
            elseif isequal(ClinicalData.Histology{s},'Non-seminoma - pure')
                patch(x1,y1,[255,51,153]./255,'edgecolor','none')
            
            elseif isequal(ClinicalData.Histology{s},'Non-seminoma - pure (teratoma)') || isequal(ClinicalData.Histology{s},'Non-seminoma - TSMT') || isequal(ClinicalData.Histology{s},'Non-seminoma - mixed (teratoma predominant)')
                patch(x1,y1,[102 204 102]./255,'edgecolor','none')
            elseif isequal(ClinicalData.Histology{s},'Other')
                patch(x1,y1,[51,102,153]./255,'edgecolor','none')
            end
        end
    
    
  
    if v==4
        if isequal(ClinicalData.Testes_Sample{s},'TGCT')
        
            patch(x1,y1,[.5 .5 .5],'edgecolor','none')
        elseif isequal(ClinicalData.Testes_Sample{s},'TGCT-MET')
            patch(x1,y1,[0 0 0],'edgecolor','none')
        elseif isequal(ClinicalData.Testes_Sample{s},'GCT-Mediastinal')
                patch(x1,y1,[1 0 0],'edgecolor','none')
        end
    end
    end

end
Labels = {'OS','Resistance','Histology','Location','KRAS*','RPL5*','TP53','12p','KRAS amp'};
for v=1:length(Events)+extra_tracks
    j=NV-v+extra_tracks;
    text (0,j-.15,Labels{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end
M_ABS.ccf_hat=str2double(M_ABS.ccf_hat);

figure()
M_ABS=reorder_struct(M_ABS,~ismember(M_ABS.i_failure_reasons,'fstar_tumor_lod'));
for i=1:slength(ClinicalData)
clonal_subclonal_rates(i,1)=sum(M_ABS.ccf_hat(ismember(M_ABS.sample,ClinicalData.sample(i)))>.9)/28.834;
clonal_subclonal_rates(i,2)=sum(M_ABS.ccf_hat(ismember(M_ABS.sample,ClinicalData.sample(i)))<=.9)/28.834;
end
bar(clonal_subclonal_rates,'stack')
