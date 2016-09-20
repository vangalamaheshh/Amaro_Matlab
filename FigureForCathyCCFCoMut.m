M=load_struct('/Volumes/xchip_cga_home/amaro/CLL/FinalMafFreeze/final_538_MAF_v1.txt');
sig_genes=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/1_09_2015_PoNCut4_WithSalvage/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt');
sig_genes=reorder_struct(sig_genes,str2double(sig_genes.q)<=.1);
sig_genes=reorder_struct(sig_genes,~ismember(sig_genes.gene,{'ADAM30';'NHS';'TTK'})); %genes which are not expressed / did not validate
CCF_color_grade=flipud(hot(100));
code=get_coding_class_muts();
M=reorder_struct(M,ismember(M.Variant_Classification,code));
N=length(unique(M.sample));
load('EventOrderCLLSF3B1.mat');

SEG=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/SCNA_CCF_recall_5.17/CCF_aggregate_called.seg');
SEG=reorder_struct(SEG,ismember(SEG.Chromosome,{'11','12','17','13'}));
SEG.CCF_hat=str2double(SEG.CCF_hat);
%%%%

%clear sig_genes
%sig_genes.gene={'SF3B1';'IKZF3';'RPS15'};
previous_table=load_table('/Users/amaro/Downloads/matrix_clonal_538_v3.txt');
CNV_IGHV=load_table('~/Documents/CLL_Matrix_SF3B1.txt');
CNV_IGHV=rmfield(CNV_IGHV,{'header','headline'});
CNV_IGHV=reorder_struct(CNV_IGHV,CNV_IGHV.SF3B1==1);

for i=1:slength(CNV_IGHV)
    if CNV_IGHV.del11q(i)>0
        CNV_IGHV.ccf_del11q(i,1)=max(SEG.CCF_hat(ismember(SEG.sample,CNV_IGHV.seg_id{i})&ismember(SEG.Chromosome,'11')));
    else
        CNV_IGHV.ccf_del11q(i,1)=0;
    end
    if CNV_IGHV.tri_12(i)>0
       CNV_IGHV.ccf_tri_12(i,1)=max(SEG.CCF_hat(ismember(SEG.sample,CNV_IGHV.seg_id{i})&ismember(SEG.Chromosome,'12')));

    else
CNV_IGHV.ccf_tri_12(i,1)=0;
    end
    if CNV_IGHV.del13q(i)>0
    
      CNV_IGHV.ccf_del13q(i,1)=max(SEG.CCF_hat(ismember(SEG.sample,CNV_IGHV.seg_id{i})&ismember(SEG.Chromosome,'13')));
else
CNV_IGHV.ccf_del13q(i,1)=0;
    end
    if CNV_IGHV.del17p(i)>0
      CNV_IGHV.ccf_del17p(i,1)=max(SEG.CCF_hat(ismember(SEG.sample,CNV_IGHV.seg_id{i})&ismember(SEG.Chromosome,'17')));

    else
 CNV_IGHV.ccf_del17p(i,1)=0;
    end
end


M=reorder_struct(M,ismember(M.sample,CNV_IGHV.samples(ismember(CNV_IGHV.SF3B1,1))));
CNV_IGHV=reorder_struct(CNV_IGHV,ismember(CNV_IGHV.SM_ids,M.sample(ismember(M.Hugo_Symbol,'SF3B1'))));
sig_genes=reorder_struct(sig_genes,ismember(sig_genes.gene,M.Hugo_Symbol));
code=get_coding_class_muts;
M=reorder_struct(M,ismember(M.Variant_Classification,code));
%%%%%%%
NV=slength(sig_genes);
for i=1:slength(sig_genes)
    NumEvents(i)=sum(ismember(M.Hugo_Symbol,sig_genes.gene{i}));
end


figure('position',[100 100 1000 1200])
hold on
subplot('position',[0.89 0.05 .1 0.8])
h=barh(fliplr(NumEvents./N),'Facecolor',[241/255,102/255,104/255],'edgecolor',[0 0 0],'Linestyle','none');
ylim([.55 3.45])
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

M=reorder_struct(M,ismember(M.Hugo_Symbol,sig_genes.gene));
N=length(unique(M.sample));



IDS=unique(M.sample);

% for i=1:slength(sig_genes)
%     for j=1:length(IDS)
%     CNV_IGHV.IDS{j,1}=IDS{j};
%     CNV_IGHV.(sig_genes.gene{i})(j,1)=ismember(IDS{j},M.coded_id(ismember(M.Hugo_Symbol,sig_genes.gene{i})));
%
%     end
% end
ff=fieldnames(CNV_IGHV);

NV=slength(sig_genes)+1;
CNV_IGHV=sort_struct(CNV_IGHV,events,repmat(-1,length(events),1));


figure('position',[100 100 1000 1200])
hold on

set(gca,'visible','off');
%setting dims of boxs
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;


fieldstokeep={ff{ismember(ff,sig_genes.gene)|ismember(ff,{'del11q','ccf_del11q','samples','tri_12','del17p','del13q','ccf_tri_12','ccf_del17p','ccf_del13q'})}};
matrix=keep_fields(CNV_IGHV,fieldstokeep);
for s=1:slength(matrix)
    matrix.SF3B1(s,1)=max(str2double(M.ccf_median_hack(ismember(M.sample,matrix.samples{s})&ismember(M.Hugo_Symbol,'SF3B1'))));
end
events{3}='ccf_del11q';events{2}='ccf_del13q';events{7}='ccf_del17p';events{9}='ccf_tri_12';

matrix=sort_struct(matrix,{events{2:end}},repmat(-1,length({events{2:end}}),1));
events{3}='del11q';events{2}='del13q';events{7}='del17p';events{9}='tri_12';
ff2=fieldnames(matrix);
NV=length(events);
ylim([-.5,44])
for v=1:NV
    for s=1:N
        %these correspond to the shapes of the panels in this case rectangles
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
        
        patch(x1,y1,[1 1 1],'edgecolor',[0 0 0])
        
        
    end
end
sig_genes.gene{41}='11q';sig_genes.gene{42}='tri12';sig_genes.gene{43}='del13q';sig_genes.gene{44}='del17p';
for v=1:NV
    event=events{v};
    for s=1:N
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        j=NV-v;
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
        c_sample=matrix.samples{s};
        if ~isempty(find(ismember(M.sample,matrix.samples{s})&ismember(M.Hugo_Symbol,events{v}))) && ~isequal(event,'del11q')
            patch(x1,y1,CCF_color_grade(round(100*max(str2double(M.ccf_median_hack(ismember(M.sample,matrix.samples{s})&ismember(M.Hugo_Symbol,events{v}),1)))),:),'edgecolor',[0 0 0])
        elseif matrix.(event)(s)==1 && isequal(event,'del11q')
            patch(x1,y1,CCF_color_grade(round(100*matrix.ccf_del11q(s)),:),'edgecolor',[0 0 0])
elseif matrix.(event)(s)>0 && isequal(event,'del13q')
    patch(x1,y1,CCF_color_grade(round(100*matrix.ccf_del13q(s)),:),'edgecolor',[0 0 0])
elseif matrix.(event)(s)>0 && isequal(event,'del17p')
           patch(x1,y1,CCF_color_grade(round(100*matrix.ccf_del17p(s)),:),'edgecolor',[0 0 0]) 
elseif matrix.(event)(s)>0 && isequal(event,'tri_12')
patch(x1,y1,CCF_color_grade(round(100*matrix.ccf_tri_12(s)),:),'edgecolor',[0 0 0]) 
        end
    end
end

for v=1:length(events)
    j=NV-v-.75;
    text (0,j,events{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end
