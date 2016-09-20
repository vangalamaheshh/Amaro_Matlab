% CLL landscape figure
% featuring mutations,CNVs,IGHV,prior treatment,Relapse_sample_or_not
%clear
JustStilgenbauer=0;


M=load_struct('/Volumes/xchip_cga_home/amaro/CLL/FinalMafFreeze/final_538_MAF_v1.txt');
CNV_IGHV=load_table('/Users/amaro/Downloads/matrix_clonal_538_v3.txt');

sample_map=load_struct('/Users/amaro/Downloads/map_160_dfci_sm_ids.txt');
[i m]=ismember(sample_map.old,CNV_IGHV.samples);
CNV_IGHV.samples(m(m>0))=sample_map.new(i);
CNV_IGHV=rmfield(CNV_IGHV,{'header','headline'});
sample_map=load_struct('/Volumes/xchip_cga_home/amaro/CLL/SampleTableICGCStilgenbauerDFCI.txt');
pair_map=load_struct('/Volumes/xchip_cga_home/amaro/CLL/PairTableICGCStilgenbauerDFCI.txt');
external_id=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/external_id_captureStil');
for i=1:slength(external_id)
external_id.ind{i,1}=external_id.sample_id{i}(1:13);
end
[nn j]=count(external_id.ind);
CNV_IGHV.Relapse=ismember(CNV_IGHV.samples,j(nn>2));
CNV_IGHV.DataSet=ones(slength(CNV_IGHV),1)*4;
CNV_IGHV.DataSet(~cellfun(@isempty,strfind(CNV_IGHV.samples,'GCLL')),1)=3;
CNV_IGHV.DataSet(~cellfun(@isempty,strfind(CNV_IGHV.samples,'ICGC')),1)=2;



for i=1:slength(pair_map)
pair_map.indi{i,1}=pair_map.pair_id{i}(1:13);
end

[i m]=ismember(pair_map.indi,CNV_IGHV.samples);
CNV_IGHV.samples(m(m>0))=pair_map.case_sample(i);
CNV_IGHV=reorder_struct(CNV_IGHV,~ismember(CNV_IGHV.samples,{'CLL-GCLL-0205-Tumor-SM-4DP8G';'CLL-CW145-Tumor-SM-24AC1'})); %blacklisted for <5% purity

if JustStilgenbauer==1
    CNV_IGHV=reorder_struct(CNV_IGHV,CNV_IGHV.DataSet==3);
    pair_map=reorder_struct(pair_map,ismember(pair_map.case_sample,CNV_IGHV.samples));
    M=reorder_struct(M,ismember(M.sample,CNV_IGHV.samples));

end



sig_genes=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/1_09_2015_PoNCut4_WithSalvage/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt');
sig_genes.q=str2double(sig_genes.q);

sig_genes=reorder_struct(sig_genes,sig_genes.q<=.1);
sig_genes=reorder_struct(sig_genes,~ismember(sig_genes.gene,{'ADAM30';'NHS';'TTK'})); %genes which are not expressed / did not validate 
%blacklist_sites=load_struct('/Users/amaro/Downloads/salvage_IGV_no_comment.txt');
%blacklist.key=strcat(blacklist_sites.start,blacklist_sites.samples);
M.key=strcat(M.Start_position,M.sample);
%M=reorder_struct(M,~ismember(M.key,blacklist.key));

f=fieldnames(CNV_IGHV);
CNevents={f{2:end-5}};
CNevents={CNevents{45:55}};
for i=1:length(CNevents)
    CNV_IGHV.(CNevents{i})(CNV_IGHV.(CNevents{i})>0)=1;
    Ne(i,1)=sum((CNV_IGHV.(CNevents{i})));
end
Msig=reorder_struct(M,ismember(M.Hugo_Symbol,sig_genes.gene));
[j n]=count(Msig.Hugo_Symbol); %sorted alphabetically number of events per gene
Events_Table.Events=[CNevents';unique(sig_genes.gene)]; %just to sort it alphabetically
Events_Table.NumEvents=[Ne;j];
Events_Table=sort_struct(Events_Table,'NumEvents',-1);
for i=1:slength(sig_genes)
    for j=1:slength(CNV_IGHV)
    
    CNV_IGHV.(sig_genes.gene{i})(j,1)=ismember(CNV_IGHV.samples{j},Msig.sample(ismember(Msig.Hugo_Symbol,sig_genes.gene{i})));
    
    end
end

CNV_IGHV=sort_struct(CNV_IGHV,Events_Table.Events,repmat(-1,slength(Events_Table),1));


N=slength(CNV_IGHV);
NV=slength(Events_Table);
figure('position',[100 100 1000 1200])
hold on
subplot('position',[0.89 0.05 .1 0.8])
h=barh(flipud(Events_Table.NumEvents)./slength(CNV_IGHV),'Facecolor',[241/255,102/255,104/255],'edgecolor',[0 0 0],'Linestyle','none');
ylim([-2.33 slength(Events_Table)])
baseline1 = get(h,'BaseLine');
set(baseline1,'LineStyle','none');
set(gca,'YColor','white','YTick',[0 .15 .3 .45])
hold on
subplot('position',[0.05 0.05 .82 0.82])
%axis([-5 N+2 -2 NV+1]);
set(gca,'visible','off');
%setting dims of boxs 
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;

for v=1:NV
    for s=1:N
     %these correspond to the shapes of the panels in this case rectangles   
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-v;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.95 .95 .95],'edgecolor',[.95 .95 .95])

        
    end
end

amp=[1 0 0];
del=[0 0 1];
missense=[179/255 222/255 105/255];
frameshift=[102/255 0/255 204/255];
nonsense=[153/255 0 13/255];
splice=[1 1 51/255];
silent=[237/255 248/255 251/255];
ICGC=[228/255 26/255 28/255];
DFCI=[55/255 126/255 184/255];
Stilg=[77/255 175/255 74/255];
inframe=[1 51/255 1];
colors=[236,250,80;80,236,250;250,5,128;5,250,127];

for v=1:NV
    event=Events_Table.Events{v};
    for s=1:N
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
        c_sample=CNV_IGHV.samples{s};
       if CNV_IGHV.(event)(s)==1
        if ismember(event,Msig.Hugo_Symbol)
            %the event is a mutation
            type=Msig.Variant_Classification(ismember(Msig.Hugo_Symbol,event)&ismember(Msig.sample,c_sample));
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
        elseif strfind(event,'del')
            %the event is a del
            patch(x1,y1,del,'edgecolor',del)
        else
            %the event is an amp
            patch(x1,y1,amp,'edgecolor',amp)
        end
       end
        
    end
end


for v=1:length(Events_Table.Events)
    j=NV-v;
    text (0,j-.5,Events_Table.Events{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end

xlim([-20 N+2])
ylim([-5 NV+2])
extrack={'prior_trt';'IGHV';'Relapse';'DataSet'};
%CNV_IGHV.IGHV(CNV_IGHV.IGHV>0)=1;
for tracks=1:length(extrack)
    for s=1:N
    %adding additional tracks
     j=0-tracks;
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
     y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
     patch(x1,y1,[.95 .95 .95],'edgecolor',[.95 .95 .95])
     if CNV_IGHV.(extrack{tracks})(s)==1
         patch(x1,y1,[0 0 0],'edgecolor',[0 0 0])
     end
     if CNV_IGHV.(extrack{tracks})(s)==4
         patch(x1,y1,[DFCI],'edgecolor',[DFCI])
     end
     if CNV_IGHV.(extrack{tracks})(s)==3
         patch(x1,y1,[Stilg],'edgecolor',[Stilg])
     end
     if CNV_IGHV.(extrack{tracks})(s)==2
         patch(x1,y1,[ICGC],'edgecolor',[ICGC])
     end
    
    end
    text (0,j-.5,extrack{tracks},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end
%tightfig
