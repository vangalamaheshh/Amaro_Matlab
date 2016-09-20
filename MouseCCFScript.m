


%P=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.29Apr2013_allelic_CAPSEG_04.30.PP-calls_tab.29May2013.txt')%Ewings PP call chart

%gender_info_file=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.sif.ABSOLUTE.txt');

P=load_struct('~/MouseInfo.txt')

%SETD=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/DiagnosticEwingSet.txt') %replace with Ewings Diagnostic pair info
%SETT=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/MetRecurEwingsSet.txt') %replace with relapse pair info
SEGS=load_struct('~/Projects/Mouse_SCLC/BarcodePassSegs/MouseBarcodeSegs.seg');
%SEGS=load_struct('/xchip/cga/gdac-prod/cga/jobResults/CapSegModule/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN/2735575/0.CapSegModule.Finished/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.segments.tsv');
%XD=load_table(''); %replace with Ewings Diagnostic MAF
%XT=load_table(''); %replace with relapse MAF
%X0=load_struct('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN/2830114/0.MutSigRun.Finished/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.final_analysis_set.maf');

 X0=load_struct('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/QC_pass_BarcodeSet/3952388/0.MutSigRun.Finished/QC_pass_BarcodeSet.final_analysis_set.maf');
 
 for i=1:slength(SEGS)
SEGS.individual_id{i,1}=strrep(SEGS.Sample{i},'.','-');
 end
 individuals=unique(SEGS.individual_id);
 individuals_mut=individuals;
 individuals_mut{25}='MOUSE_SCLC-AD3588-TM-NB-SM-1JIOR-SM-1JINT';
 individuals_mut{26}='MOUSE_SCLC-AD3588-TM-NB-SM-2YUPG-SM-1JINT';
 individuals_mut{27}='MOUSE_SCLC-AD3588-TP-NB-SM-1JIOO-SM-1JINT';
 individuals_mut{28}='MOUSE_SCLC-AD3588-TP-NB-SM-1JIOP-SM-1JINT';
 individuals_mut{29}='MOUSE_SCLC-AD3588-TP-NB-SM-1JIOQ-SM-1JINT';
 X0_rejects=reorder_struct(X0,~ismember(X0.patient,individuals_mut));
 X0=reorder_struct(X0,ismember(X0.patient,individuals_mut));
 
for i=1:length(individuals)
MUTS=ismember(X0.patient,individuals_mut{i});
MUTS=find(MUTS==1);

SEG=reorder_struct(SEGS,ismember(SEGS.individual_id,individuals{i}));
SEG.End=cellfun(@str2num,SEG.End);
SEG.Start=cellfun(@str2num,SEG.Start);

    for j=1:length(MUTS)
    pos=X0.Start_position(MUTS(j));
    pos=cell2mat(pos);
    chr=X0.Chromosome{MUTS(j)};
    if isequal(chr,'X') || isequal(chr,'Y') || isequal(chr,'M')
        X0.segvalue{MUTS(j),1}='0';
       
    else
    locs=find(ismember(SEG.Chromosome,sprintf('chr%s',chr))==1);
    ll=find(SEG.End(locs)>pos(1),1,'first');
    X0.segvalue{MUTS(j),1}=SEG.Segment_Mean{ll};
 
    
    end
    end

    
end

X0.segvalue=cellfun(@str2num,X0.segvalue);


CCF_range=0:0.01:1;
%X0.SET=repmat({''},length(X0.patient),1);
%X0.SET=X0.SET';
%kD=find(ismember(X0.patient,SETD.pair_id));
%X0.SET(kD)={'Primary'};
%kT=find(ismember(X0.patient,SETT.pair_id));
%X0.SET(kT)={'Recurr'};
%X0.SET=X0.SET';
%z=repmat('|',length(X0.Chromosome),1);
%X0.id=regexprep(cellstr([char(X0.Chromosome) z char(X0.Start_position) z char(X0.Tumor_Sample_Barcode)]),' ','');

X0.purity=NaN*zeros(size(X0.patient));

%X0.homozygous_ix=NaN*zeros(size(X0.id));
X0.cancer_cell_frac=NaN*zeros(size(X0.patient));
k=find(isnan(X0.cancer_cell_frac));

for i=1:slength(X0)
strs=regexp(X0.patient{i},'-','split');
X0.pat{i,1}=strs{6};
end
for i=1:slength(P)
P.pat{i,1}=P.Tumor{i}(end-4:end);
end



[i m]=ismember(X0.pat,P.pat);
P0=trimStruct(P,m(m>0));
P0.ploidy=ones(size(P0.pat)).*2;
X0.gender=P0.gender;
P0.Purity=cellfun(@str2num,P0.Purity,'UniformOutput', false);
P0.Purity=cell2mat(P0.Purity);
%%% correcting for ploidy in purity and %

var=(P0.Purity./100)./(1-(P0.Purity./100));
var2=(2./P0.ploidy);
P0.PurityNew=(var.*var2)./(var.*var2+1)
X0.purity(k)=P0.PurityNew(k);
X0.ploidy=P0.ploidy;
X0=reorder_struct(X0,X0.purity>0);
k=find(isnan(X0.cancer_cell_frac));

for i=1:length(X0.segvalue)
    
X0.capseg_val(i,1)=(2^X0.segvalue(i)*2)/X0.purity(i);

end



X0.t_ref_count=cellfun(@str2num,X0.t_ref_count);
X0.t_alt_count=cellfun(@str2num,X0.t_alt_count);
X0.i_tumor_f=cellfun(@str2num,X0.i_tumor_f);


X0.O_tumor=X0.t_alt_count./(X0.t_alt_count+X0.t_ref_count);


X0.CN=X0.capseg_val.*(X0.ploidy./2);
%Xs=ismember(X0.Chromosome,'X');
%Males=ismember(X0.gender,'Male');
%MX=find((Xs+Males)==2);
%X0.CN(MX)=X0.CN(MX)./2;

for i=1:length(X0.t_alt_count)
    
   if isequal(X0.Chromosome{i},'X') && isequal(X0.gender{i},'M')
           
           CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
           [val1,tmo1]=max(CCFPdist1);
           X0.CCF(i,1)=CCF_range(tmo1);
[X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
   else
    
   CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
   [val1,tmo1]=max(CCFPdist1);
   CCFPdist2=betapdf(2.*(CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
   [val2,tmo2]=max(CCFPdist2);
   if val1>val2 || X0.CN(i)<1.5
[X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
 X0.CCF(i,1)=CCF_range(tmo1);
   else
[X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
 X0.CCF(i,1)=CCF_range(tmo1);
   end
   end
   
end

for i=1:length(X0.pci)
    X0.total_count(i,1)=X0.t_alt_count(i)+X0.t_ref_count(i);
    X0.af_CI95_low(i,1) = X0.t_alt_count(i)-(X0.pci{i}(1)*X0.total_count(i,1));
    X0.af_CI95_high(i,1) = X0.pci{i}(2)*X0.total_count(i,1)+X0.t_alt_count(i); 
    
    if X0.af_CI95_high(i)>X0.total_count(i)
        X0.af_CI95_high(i)=X0.total_count(i);
    end
    
    CCFdist_low=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.af_CI95_low(i)+1,X0.t_ref_count(i)+1);
    CCFdist_high=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.af_CI95_high(i)+1,X0.t_ref_count(i)+1);
    [~,tmo_low]=max(CCFdist_low);
    [~,tmo_high]=max(CCFdist_high);
    X0.high_CCF_CI(i,1)=CCF_range(tmo_high);
    X0.low_CCF_CI(i,1)=CCF_range(tmo_low);
end

pats=unique(X0.patient);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% optimizing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ploidy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for p=1:length(unique(X0.patient))
%     pats{p}
%  [nelements,xcenters]=hist(X0.CCF(ismember(X0.patient,pats{p})));
%  
%  while max(xcenters)<.95
% 
%      X0.ploidy(ismember(X0.patient,pats{p}))=X0.ploidy(ismember(X0.patient,pats{p}))+.01;
%      X0.CN=X0.capseg_val.*(X0.ploidy./2);
% 
%      median(X0.ploidy(ismember(X0.patient,pats{p})))
%      max(xcenters)
%      median(X0.purity(ismember(X0.patient,pats{p})))
%     var=(X0.purity)./(1-(X0.purity));
%     var2=(2./X0.ploidy);
%     X0.purity=(var.*var2)./(var.*var2+1);
%     
%     for i=1:length(X0.t_alt_count)
%     
%    if isequal(X0.Chromosome{i},'X') && isequal(X0.gender{i},'M')
%            
%            CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%            [val1,tmo1]=max(CCFPdist1);
%            X0.CCF(i,1)=CCF_range(tmo1);
%     [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%    else
%     
%    CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%    [val1,tmo1]=max(CCFPdist1);
%    CCFPdist2=betapdf(2.*(CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%    [val2,tmo2]=max(CCFPdist2);
%    if val1>val2 || X0.CN(i)<1.5
%     [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%     X0.CCF(i,1)=CCF_range(tmo1);
%    else
%     [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%     X0.CCF(i,1)=CCF_range(tmo1);
%    end
%    end
%    
%    end
% 
%     
%      
%    for i=1:length(X0.pci)
%     X0.total_count(i,1)=X0.t_alt_count(i)+X0.t_ref_count(i);
%     X0.af_CI95_low(i,1) = X0.t_alt_count(i)-(X0.pci{i}(1)*X0.total_count(i,1));
%     X0.af_CI95_high(i,1) = X0.pci{i}(2)*X0.total_count(i,1)+X0.t_alt_count(i); 
%     
%     if X0.af_CI95_high(i)>X0.total_count(i)
%         X0.af_CI95_high(i)=X0.total_count(i);
%     end
%     
%     CCFdist_low=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.af_CI95_low(i)+1,X0.t_ref_count(i)+1);
%     CCFdist_high=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.af_CI95_high(i)+1,X0.t_ref_count(i)+1);
%     [~,tmo_low]=max(CCFdist_low);
%     [~,tmo_high]=max(CCFdist_high);
%     X0.high_CCF_CI(i,1)=CCF_range(tmo_high);
%     X0.low_CCF_CI(i,1)=CCF_range(tmo_low);
%    end
% 
%      
%      [nelements,xcenters]=hist(X0.CCF(ismember(X0.patient,pats{p})));
%  
%  end
% end

X0.low_CCF_CI(X0.low_CCF_CI==1)=.9;
X0.subclonal_ix(k,1)=X0.CCF(k)<0.9;
save_struct(X0,'/xchip/cga_home/amaro/Mouse_SCLC/CCF_Maf_MatlabBarcodeSet.txt')

%%%%%%%%%%%%%% Clustering

X=X0;
X.cancer_cell_frac = X.CCF;
z=repmat('|',length(X0.Chromosome),1);
X.site=regexprep(cellstr([char(X0.Chromosome) z char(X0.Start_position) z char(X0.Tumor_Seq_Allele1)]),' ','');
X.idp=regexprep(strcat(X.Hugo_Symbol,':',X.Protein_Change,':',X.Chromosome,':',cellstr(char(X.Start_position))),' ','');
X.id1=regexprep(strcat(X.Chromosome,':',cellstr(char(X.Start_position)),'_',X.Tumor_Sample_Barcode),' ','');

%%%%%%%%% Number of clones for each sample table
pats=unique(X.patient);
for i=1:length(pats)
    
   
    pat_maf=reorder_struct(X,ismember(X.patient,pats{i}));
    distance=(pat_maf.high_CCF_CI-pat_maf.low_CCF_CI)/2;
    pat_maf.CloneAssignment=cluster1D(pat_maf.CCF,distance)';
    table.patient{i,1}=pats{i,1};
    table.NumberOfSubclones(i,1)=max(pat_maf.CloneAssignment);
    
end

save_struct(table,'~/Projects/Mouse_SCLC/NumberOfSubclones.txt');

pair_mapping_file=load_struct('~/Projects/Mouse_SCLC/pair_map_CCF.txt');

CMR=1;
for zz=1:length(pair_mapping_file.pair_id_p)
    clear vals1 vals2 vals3 color_vector1 color_vector2 color_vector3
    pri_ind= ismember(X.patient, pair_mapping_file.pair_id_p{zz});
    met_ind = ismember(X.patient, pair_mapping_file.pair_id_m{zz});
    x1 = reorder_struct(X, pri_ind);
    x2 = reorder_struct(X, met_ind);
clear xy xyclust yc xc rx ry counter 

counter=1;
 for i=1:slength(x2)
 if sum(ismember(x1.idp,x2.idp{i}))
 xy(i,1)=x1.CCF(find(ismember(x1.idp,x2.idp{i})==1));
 xy(i,2)=x2.CCF(i);
 xy(i,3)=i;
 xy(i,4)=x1.high_CCF_CI(find(ismember(x1.idp,x2.idp{i})==1));
 xy(i,5)=x1.low_CCF_CI(find(ismember(x1.idp,x2.idp{i})==1));
 xy(i,6)=x2.high_CCF_CI(i);
 xy(i,7)=x2.low_CCF_CI(i);
 xyclust(counter,1)=x1.CCF(find(ismember(x1.idp,x2.idp{i})==1));
 xyclust(counter,2)=x2.CCF(i);
 xyclust(counter,4)=x2.high_CCF_CI(i)-x2.low_CCF_CI(i);
 xyclust(counter,3)=x1.high_CCF_CI(find(ismember(x1.idp,x2.idp{i})==1))-x1.low_CCF_CI(find(ismember(x1.idp,x2.idp{i})==1));
 counter=counter+1;
 else
 xy(i,1)=0;
 xy(i,2)=x2.CCF(i);
 xy(i,3)=0;
 
 
 xy(i,4)=0.008;
 xy(i,5)=0;
 xy(i,6)=x2.high_CCF_CI(i);
 xy(i,7)=x2.low_CCF_CI(i);
 
 
 
 end
 end

j=length(xy);
counter=1;
for i=1:slength(x1)
    if ~sum(ismember(x2.idp,x1.idp{i}))
        xy(counter+j,1)=x1.CCF(i);
        xy(counter+j,2)=0;
        xy(counter+j,3)=0;
        xy(counter+j,4)=x1.high_CCF_CI(i);
        xy(counter+j,5)=x1.low_CCF_CI(i);
        xy(counter+j,6)=0.008;
        xy(counter+j,7)=0;
        counter=counter+1;

    end
end
min_distance=.6;
G=clusterNN(xyclust(:,1:2),xyclust(:,3:4)*min_distance);
xy(find(xy(:,3)~=0),3)=G;
next_cluster=max(G)+1;
for i=1:length(xy)
   if xy(i,3)==0 && xy(i,2)==0
       xy(i,3)=next_cluster;
       
   end
   if xy(i,3)==0 && xy(i,1)==0
       xy(i,3)=next_cluster+1;
   end

end




% t=0:.01:2*pi;
for i=1:max(xy(:,3))
%     if i<next_cluster
     pts=xy(find(xy(:,3)==i),1:7);
     %1/sqrt(sum(1/CI^2))
        CIx=pts(:,4)-pts(:,5);
        CIy=pts(:,6)-pts(:,7);
        CIx=CIx.*CIx;
        CIy=CIy.*CIy;
        
        
        rx(i,1)=1/sqrt(sum(1/CIx));
        ry(i,1)=1/sqrt(sum(1/CIy));
    xc(i,1)=median(pts(:,1));
     yc(i,1)=median(pts(:,2));
%     plot(r*sin(t)+yc,r*cos(t)+xc);
%     hold on
%     else
%     pts=xy(find(xy(:,3)==i),1:2);
%     r=max(max(dist(pts(:,1),pts(:,2),'euclidean'))/4)
%     xc=mean(pts(:,1));
%     yc=mean(pts(:,2));
%     end
end

for i=1:length(xc)
CM.CCF(CMR,1)=xc(i,1);
CM.patient{CMR,1}=x1.patient{1};
CM.high_CCF_CI(CMR,1)=CM.CCF(CMR,1)+rx(i,1)/2;
CM.low_CCF_CI(CMR,1)=CM.CCF(CMR,1)-rx(i,1)/2;


CMR=1+CMR;
end

for i=1:length(yc)
CM.CCF(CMR,1)=yc(i,1);
CM.patient{CMR,1}=x2.patient{1};
CM.high_CCF_CI(CMR,1)=CM.CCF(CMR,1)+ry(i,1)/2;
CM.low_CCF_CI(CMR,1)=CM.CCF(CMR,1)-ry(i,1)/2;
CMR=CMR+1;
end

end
% 
save_struct(CM,'~/Projects/Mouse_SCLC/ClustersMAF.txt');


% [i m]=ismember(x1.idp,x2.idp);
% x12=reorder_struct(x1,find(i));
% x21=reorder_struct(x2,m(m>0));
% x10=reorder_struct(x1,find(~i));
% [i m]=ismember(x2.idp,x1.idp); 
% x20=reorder_struct(x2,find(~i));
% 
% %k=find((x10.cancer_cell_frac>0.7)&(x10.i_tumor_f>0.7));
% %k=find((x10.cancer_cell_frac>0.7));
% %fprintf('%s\n',x10.id1{k})
% %k=find((x20.cancer_cell_frac>0.7));
% %fprintf('%s\n',x20.id1{k})
% 
% clf
% nc=slength(x12);
% %nc = x12.N
% a=0:359;
% cm=jet(nc);
% for i=1:size(x12.CCF)
% vals1{i}=sprintf('%d-%d',x12.CCF(i),x21.CCF(i));
% 
% end
% 
% if length(x12.CCF)==0
%     vals1=[];
%     color_vector1=[];
% else
% 
% [cs,strs]=count(vals1);
% for i=1:size(strs)
%    color_vector1(ismember(vals1,strs{i}))=cs(i);
% end
% end
%  %plot(x12.cancer_cell_frac,x21.cancer_cell_frac,'o');
%  %text(x12.cancer_cell_frac+0.02,x21.cancer_cell_frac+0.02,x21.Hugo_Symbol,'verticalalignment','top','fontsize',9);
% 
% 
% for i=1:nc
% 
%     x0=x12.cancer_cell_frac(i);
%     y0=x21.cancer_cell_frac(i);
% 
%     x1=cosd(a);
%     x2=0*x1;
%    
%     k=find((x1>0));
%     x2(k)=x1(k)*(x12.high_CCF_CI(i)-x0);
%     k=find((x1<=0));
%     x2(k)=x1(k)*(x0-x12.low_CCF_CI(i));
% 
%     y1=sind(a);
%     y2=0*y1;
%     k=find((y1>0));
%     y2(k)=y1(k).*(x21.high_CCF_CI(i)-y0);
%     k=find((y1<=0));
%     y2(k)=(y1(k)).*(y0-x21.low_CCF_CI(i));
%     patch(x0+x2,y0+y2,cm(i,:),'edgecolor','none');
% end
% 
% hold on;
% nc=slength(x10);
% cm=jet(nc);
% for i=1:size(x10.cancer_cell_frac)
% vals2{i}=sprintf('%d-0',x10.cancer_cell_frac(i));
% end
% [cs,strs]=count(vals2);
% 
% for i=1:size(strs)
%    color_vector2(ismember(vals2,strs{i}))=cs(i);
% end
% 
% %plot(x10.cancer_cell_frac,0*x10.cancer_cell_frac,'k','.')
% %text(x10.cancer_cell_frac+0.02,0*x10.cancer_cell_frac+0.02,x10.Hugo_Symbol,'verticalalignment','top','fontsize',9)
% 
% for i=1:nc
% 
%     x0=x10.cancer_cell_frac(i);
%     y0=0*x10.cancer_cell_frac(i);
% 
%     x1=cosd(a);
%     x2=0*x1;
%     k=find((x1>0));
%     x2(k)=x1(k).*(x10.high_CCF_CI(i)-x0);
%     k=find((x1<=0));
%     x2(k)=x1(k).*(x0-x10.low_CCF_CI(i));
% 
%     y1=sind(a);
%     y2=y1*0.01;
%     patch(x0+x2,y0+y2,cm(i,:),'edgecolor','none');
% end
% 
% nc=slength(x20);
% cm=jet(nc);
% 
% for i=1:size(x20.cancer_cell_frac)
% vals3{i}=sprintf('%d-0',x20.cancer_cell_frac(i));
% end
% [cs,strs]=count(vals3);
% 
% for i=1:size(strs)
%    color_vector3(ismember(vals3,strs{i}))=cs(i);
% end
% 
% Xscat=vertcat(x12.cancer_cell_frac,x10.cancer_cell_frac,0*x20.cancer_cell_frac);
% Yscat=vertcat(x21.cancer_cell_frac,0*x10.cancer_cell_frac,x20.cancer_cell_frac);
% Cscat=vertcat(color_vector1',color_vector2',color_vector3');
% Cs=(Cscat>1);
% 
% %plot(0*x20.cancer_cell_frac,x20.cancer_cell_frac,'k','.')
% %text(0*x20.cancer_cell_frac+0.02,x20.cancer_cell_frac+0.02,x20.Hugo_Symbol,'verticalalignment','top','fontsize',9)
% for i=1:nc
% 
%     x0=0*x20.cancer_cell_frac(i);
%     y0=x20.cancer_cell_frac(i);
% 
%     x1=cosd(a);
%     x2=0*x1;
%     x2=x1*0.01;
% 
%     y1=sind(a);
%     y2=0*y1;
%     k=find((y1>0));
%     y2(k)=y1(k).*(x20.high_CCF_CI(i)-y0);
%     k=find((y1<=0));
%     y2(k)=y1(k).*(y0-x20.low_CCF_CI(i));
%     patch(x0+x2,y0+y2,cm(i,:),'edgecolor','none');
% end
% scatter(Xscat,Yscat,60,Cs,'fill')
% colormap([0.75,0.75,0.75; 0,0,0])
% hPatch = findobj(gca,'Type','patch');
% set(hPatch,'facealpha',0.08)
% axis equal;
% axis([-.01 1.1 -.01 1])
% title(pair_mapping_file.individual_id{zz},'FontSize',20)
% xlabel(pair_mapping_file.pair_id_p{zz},'FontSize',20)
% ylabel(pair_mapping_file.pair_id_m{zz},'FontSize',20)
% f=sprintf('/xchip/cga_home/amaro/EwingsSarcoma/%s_EwingsCloud.png',pair_mapping_file.individual_id{zz});
% %print_D(f,{{'fig'},{'pdf'},{'png','-r300'}});
% 
% %saveas(gcf,f,'eps2c')
% print(gcf, '-dpng', '-r400',[f,'.png'])
% % frame=getframe(gcf);
% % if isempty(frame.colormap)
% %    imwrite(frame.cdata, f)
% % else
% %    imwrite(frame.cdata, frame.colormap, f,'Quality',100)
% % end
% 
% hold off
% % 
% % end
% end
% % 
% % %{%% MY OWN PLOT:
% % z=repmat('|',length(X0.Chromosome),1);
% % X0.site=regexprep(cellstr([char(X0.Chromosome) z num2str(X0.Start_position) z char(X0.Tumor_Seq_Allele1)]),' ','');
% % pri_ind= ismember(X0.patient, SETD.pair_id);
% % met_ind = ismember(X0.patient, SETT.pair_id);
% % pri = reorder_struct(X0, pri_ind);
% % met = reorder_struct(X0, met_ind);
% % 
% % [o met_ind xom_ind] = intersect(met.site, X0.site);
% % [o pri_ind xop_ind] = intersect(pri.site, X0.site);
% % 
% % total=struct;
% % total.site = X0.site;
% % total.Hugo_Symbol = X0.Hugo_Symbol
% % total.CCF_met = zeros(length(total.site),1);
% % total.CCF_pri = zeros(length(total.site),1);
% % total.CCF_met(xom_ind) = met.CCF(met_ind);
% % total.CCF_pri(xop_ind) = pri.CCF(pri_ind);
% % 
% % %Heirarchical clustering:
% % X = [total.CCF_pri total.CCF_met];
% % %Z = linkage(X,'single');
% % %T = cluster(Z, 'maxclust', 4);
% % 
% % seed = [0 0; 0 1; 1 0; 1 1];
% % [T cntr] = kmeans(X,4,'distance','cityblock','start',seed);
% % X0.clusters = T;
% % total.clusters = T;
% % 
% % clusts = unique(T);
% % col = {'ko','mo','bo','ro'}
% % 
% % for i=1:length(clusts)
% %     p = reorder_struct(total, total.clusters==clusts(i));
% %     plot(p.CCF_pri, p.CCF_met, col{i});
% %     hold on
% % end
% % 
% % xlabel('Primary CCF')
% % ylabel('Met CCF')
% % 
% 
% 
% %}
% 
% %X0=mergeStruct(XD,XT);
% 
% %X0.multiplicity=NaN*zeros(size(X0.id));
% %X0.subclonal_ix=NaN*zeros(size(X0.id));
% %X0.q_hat=NaN*zeros(size(X0.id));
% %X0.modal_q_s=NaN*zeros(size(X0.id));
% 
% %[i m]=ismember(X0.id,X.id);
% %k=find(i);
% %X0.subclonal_ix(k)=X.subclonal_ix(m(m>0));
% %X0.purity(k)=X.purity(m(m>0));
% %X0.q_hat(k)=X.q_hat(m(m>0));
% %X0.modal_q_s(k)=X.modal_q_s(m(m>0));
% %X0.cancer_cell_frac(k)=X.cancer_cell_frac(m(m>0));
% %X0.homozygous_ix(k)=X.homozygous_ix(m(m>0));
% %X0.multiplicity=X0.cancer_cell_frac.*X0.modal_q_s;
% 
%  %X0.AF_Tumor=X0.i_tumor_f./X0.purity;
%   
% %X0.multiplicity(k)=X0.i_tumor_f(k).*X0.CN(k)./X0.purity(k);
% %
% %X0.homozygous_ix(k)=round(X0.multiplicity(k))==X0.CN(k) ;
% %X0.modal_q_s(k)=round(X0.CN(k).*X0.purity(k))+ones(size(X0.CN(k)));
%  %X0.modal_q_s(k)=round((X0.CN(k)-ones(size(X0.CN(k)))*2)./X0.purity(k) + 2*ones(size(X0.CN(k))));
%  %X0.modal_q_s(X0.modal_q_s<1)=1;
% %X0.cancer_cell_frac=X0.multiplicity./X0.modal_q_s(k);
% %X0.cancer_cell_frac(find(X0.cancer_cell_frac>1))=1;
% %and to convert to CCF we can use something like this with modal_q_s:
% %X0.multiplicity=X0.cancer_cell_frac.*X0.modal_q_s;
