


P=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.29Apr2013_allelic_CAPSEG_04.30.PP-calls_tab.29May2013.txt')%Ewings PP call chart

gender_info_file=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.sif.ABSOLUTE.txt');


SETD=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/DiagnosticEwingSet.txt') %replace with Ewings Diagnostic pair info
SETT=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/MetRecurEwingsSet.txt') %replace with relapse pair info

SEGS=load_struct('/xchip/cga/gdac-prod/cga/jobResults/CapSegModule/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN/2735575/0.CapSegModule.Finished/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.segments.tsv');
%XD=load_table(''); %replace with Ewings Diagnostic MAF
%XT=load_table(''); %replace with relapse MAF
%X0=load_struct('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN/2830114/0.MutSigRun.Finished/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.final_analysis_set.maf');
X0=load_struct('~/Projects/EwingsSarcoma/An_SIGMA_EwingSarcoma_Pair.forcedCall.indel_pulldown_maf.MutSig.blacklist.19Nov2013.maf');
X0.Start_position=cellfun(@str2num,X0.Start_position);
X0.t_ref_count=cellfun(@str2num,X0.t_ref_count);
X0.t_alt_count=cellfun(@str2num,X0.t_alt_count);
X0.i_tumor_f=cellfun(@str2num,X0.i_tumor_f);
X0.id=regexprep(strcat(X0.Hugo_Symbol,'_',cellstr(num2str(X0.Start_position)),'_',X0.individual),' ','');

[counts names]=count(X0.id);
ids=names(find(counts==2));
for i=1:length(ids)
    k=find(ismember(X0.id,ids{i}));
if ismember(X0.Hugo_Symbol(k(1)),'MTHFD2L')
    
end
    if(length(unique(X0.Tumor_Seq_Allele2(k))))>1
        
    if X0.t_alt_count(k(1))>=X0.t_alt_count(k(2))
        X0.t_alt_count(k(2))=0;
        X0.Tumor_Seq_Allele2(k(2))=X0.Tumor_Seq_Allele2(k(1));
    else
        X0.t_alt_count(k(1))=0;
        X0.Tumor_Seq_Allele2(k(1))=X0.Tumor_Seq_Allele2(k(2));
    end
    end
end


for i=1:length(SEGS.individual_id)  
MUTS=ismember(X0.patient,SEGS.individual_id{i});
MUTS=find(MUTS==1);
SEG=load_struct(SEGS.capseg_segmentation_file{i});
SEG.End=cellfun(@str2num,SEG.End);
SEG.Start=cellfun(@str2num,SEG.Start);

    for j=1:length(MUTS)
    pos=X0.Start_position(MUTS(j));
    chr=X0.Chromosome{MUTS(j)};
    locs=find(ismember(SEG.Chromosome,chr));
    ll=find(SEG.End(locs)>pos,1,'first');
    
    X0.segvalue{MUTS(j),1}=SEG.Segment_Mean{locs(ll)};
 
    
    end

    
end

X0.segvalue=cellfun(@str2num,X0.segvalue);


CCF_range=0:0.01:1;
X0.SET=repmat({''},length(X0.patient),1);
X0.SET=X0.SET';
kD=find(ismember(X0.patient,SETD.pair_id));
X0.SET(kD)={'Primary'};
kT=find(ismember(X0.patient,SETT.pair_id));
X0.SET(kT)={'Recurr'};
X0.SET=X0.SET';
z=repmat('|',length(X0.Chromosome),1);
X0.id=regexprep(cellstr([char(X0.Chromosome) z char(X0.Start_position) z char(X0.Tumor_Sample_Barcode)]),' ','');

X0.purity=NaN*zeros(size(X0.id));

%X0.homozygous_ix=NaN*zeros(size(X0.id));
X0.cancer_cell_frac=NaN*zeros(size(X0.id));

k=find(isnan(X0.cancer_cell_frac));
X0.pat=X0.patient;
P.pat=P.sample_id;

[i m]=ismember(X0.pat,P.pat);
P0=trimStruct(P,m(m>0));
G0=trimStruct(gender_info_file,m(m>0));
X0.gender=G0.gender;
P0.purity=cellfun(@str2num,P0.purity,'UniformOutput', false);
P0.purity=cell2mat(P0.purity);
X0.purity(k)=P0.purity(k);
P0.ploidy=cellfun(@str2num,P0.ploidy,'UniformOutput', false);
P0.ploidy=cell2mat(P0.ploidy);
X0.ploidy=P0.ploidy;

% for i=1:length(X0.segvalue)
% X0.capseg_val(i,1)=(2^X0.segvalue(i)*2)/X0.purity(i);
% end

for i=1:length(X0.segvalue)
%         X0.capseg_val(i,1)=2*((1+(2^X0.segvalue(i)-1))/X0.purity(i));
X0.CN(i,1)=X0.ploidy(i).*(2.^X0.segvalue(i)-1+X0.purity(i))./X0.purity(i);
        %X0.capseg_val(i,1)=2*(1+(2^X0.segvalue(i)-1)/X0.purity(i));

end
X0.CN(X0.CN<0)=0;
X0.CN(X0.CN==0)=2;


%X0=trimStruct(X0,~X0.CN==0);
X0=rmfield_if_exist(X0,'header');
X0=rmfield_if_exist(X0,'headline');
X0=rmfield_if_exist(X0,'N');

X0.O_tumor=X0.t_alt_count./(X0.t_alt_count+X0.t_ref_count);


%X0.CN=X0.capseg_val.*(X0.ploidy./2);
%Xs=ismember(X0.Chromosome,'X');
%Males=ismember(X0.gender,'Male');
%MX=find((Xs+Males)==2);
%X0.CN(MX)=X0.CN(MX)./2;
for i=1:length(X0.t_alt_count)
    
   if isequal(X0.Chromosome{i},'X') && isequal(X0.gender{i},'Male')
           
           CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
           [val1,tmo1]=max(CCFPdist1);
           X0.CCF(i,1)=CCF_range(tmo1);
[X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);

   if X0.CCF(i,1)>=0.05
 for m=1:100
  if sum(CCFPdist1(1:m))/total>.10
      X0.low_CCF_CI(i,1)=CCF_range(m);
      break
  end
 end
 else
     X0.low_CCF_CI(i,1)=0;
 end
 
 if X0.CCF(i,1)<=.95
 for m=101:-1:1
  if sum(CCFPdist1(m:101))/total>.10
      X0.high_CCF_CI(i,1)=CCF_range(m);
      break
  end
 end
 else
     X0.high_CCF_CI(i,1)=1;
 end
 
   else
    
   CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
   [val1,tmo1]=max(CCFPdist1);
 total=sum(CCFPdist1);
 
   
   CCFPdist2=betapdf(2.*(CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
   [val2,tmo2]=max(CCFPdist2);
   if val1>val2 || X0.CN(i)<1.5
[X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
X0.pci_high(i,1)=X0.pci{i,1}(2);
X0.pci_low(i,1)=X0.pci{i,1}(1);
 X0.CCF(i,1)=CCF_range(tmo1);
   else
[X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
X0.pci_high(i,1)=X0.pci{i,1}(2);
X0.pci_low(i,1)=X0.pci{i,1}(1);
 X0.CCF(i,1)=CCF_range(tmo2);
 
 
 
 
 
 
   end
   
 if X0.CCF(i,1)>=0.07
 for m=1:100
  if sum(CCFPdist1(1:m))/total>.10
      X0.low_CCF_CI(i,1)=CCF_range(m);
      break
  end
 end
 else
     X0.low_CCF_CI(i,1)=0;
 end
 
 if X0.CCF(i,1)<=.93
 for m=101:-1:1
  if sum(CCFPdist1(m:101))/total>.10
      X0.high_CCF_CI(i,1)=CCF_range(m);
      break
  end
 end
 else
     X0.high_CCF_CI(i,1)=1;
 end
   
   end
end
%%%%%%%%%%%%%%%%%% THIS WAS WRONG
% for i=1:length(X0.pci)
%     X0.total_count(i,1)=X0.t_alt_count(i)+X0.t_ref_count(i);
%     X0.af_CI95_low(i,1) = X0.pci_low(i,1);
%     X0.af_CI95_high(i,1) = X0.pci_high(i,1);
%     
% 
%     
%     CCFdist_low=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.af_CI95_low(i)+1,X0.t_ref_count(i)+1);
%     CCFdist_high=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.af_CI95_high(i)+1,X0.t_ref_count(i)+1);
%     
%     [~,tmo_low]=max(CCFdist_low);
%     [~,tmo_high]=max(CCFdist_high);
%     X0.high_CCF_CI(i,1)=CCF_range(tmo_high);
%     X0.low_CCF_CI(i,1)=CCF_range(tmo_low);
% end




X0.subclonal_ix(k,1)=X0.CCF(k)<0.9;


save_struct(X0,'/xchip/cga_home/amaro/EwingsSarcoma/CCF_MAF_Ewings8.15.redone.txt')



X=X0;
X.cancer_cell_frac = X.CCF;
z=repmat('|',length(X0.Chromosome),1);
X.site=regexprep(cellstr([char(X0.Chromosome) z char(X0.Start_position) z char(X0.Tumor_Seq_Allele1)]),' ','');
X.idp=regexprep(strcat(X.Hugo_Symbol,':',X.Protein_Change,':',X.Chromosome,':',cellstr(char(X.Start_position))),' ','');
X.id1=regexprep(strcat(X.Chromosome,':',cellstr(char(X.Start_position)),'_',X.Tumor_Sample_Barcode),' ','');

pair_mapping_file=load_struct('/xchip/cga_home/amaro/EwingsSarcoma/pairmappingfileEwing.txt');


for zz=1:length(pair_mapping_file.pair_id_p)
    clear vals1 vals2 vals3 color_vector1 color_vector2 color_vector3
pri_ind= ismember(X.patient, pair_mapping_file.pair_id_p{zz});
met_ind = ismember(X.patient, pair_mapping_file.pair_id_m{zz});
x1 = reorder_struct(X, pri_ind);
x2 = reorder_struct(X, met_ind);

[i m]=ismember(x1.idp,x2.idp);
x12=reorder_struct(x1,find(i));
x21=reorder_struct(x2,m(m>0));
x10=reorder_struct(x1,find(~i));
[i m]=ismember(x2.idp,x1.idp); 
x20=reorder_struct(x2,find(~i));

%k=find((x10.cancer_cell_frac>0.7)&(x10.i_tumor_f>0.7));
%k=find((x10.cancer_cell_frac>0.7));
%fprintf('%s\n',x10.id1{k})
%k=find((x20.cancer_cell_frac>0.7));
%fprintf('%s\n',x20.id1{k})

clf
nc=slength(x12);
%nc = x12.N
a=0:359;
cm=jet(nc);
for i=1:size(x12.CCF)
vals1{i}=sprintf('%d-%d',x12.CCF(i),x21.CCF(i));

end

if length(x12.CCF)==0
    vals1=[];
    color_vector1=[];
else

[cs,strs]=count(vals1);
for i=1:size(strs)
   color_vector1(ismember(vals1,strs{i}))=cs(i);
end
end
 %plot(x12.cancer_cell_frac,x21.cancer_cell_frac,'o');
 %text(x12.cancer_cell_frac+0.02,x21.cancer_cell_frac+0.02,x21.Hugo_Symbol,'verticalalignment','top','fontsize',9);


for i=1:nc

    x0=x12.cancer_cell_frac(i);
    y0=x21.cancer_cell_frac(i);

    x1=cosd(a);
    x2=0*x1;
   
    k=find((x1>0));
    x2(k)=x1(k)*(x12.high_CCF_CI(i)-x0);
    k=find((x1<=0));
    x2(k)=x1(k)*(x0-x12.low_CCF_CI(i));

    y1=sind(a);
    y2=0*y1;
    k=find((y1>0));
    y2(k)=y1(k).*(x21.high_CCF_CI(i)-y0);
    k=find((y1<=0));
    y2(k)=(y1(k)).*(y0-x21.low_CCF_CI(i));
    patch(x0+x2,y0+y2,cm(i,:),'edgecolor','none');
end

hold on;
nc=slength(x10);
cm=jet(nc);
for i=1:size(x10.cancer_cell_frac)
vals2{i}=sprintf('%d-0',x10.cancer_cell_frac(i));
end
[cs,strs]=count(vals2);

for i=1:size(strs)
   color_vector2(ismember(vals2,strs{i}))=cs(i);
end

%plot(x10.cancer_cell_frac,0*x10.cancer_cell_frac,'k','.')
%text(x10.cancer_cell_frac+0.02,0*x10.cancer_cell_frac+0.02,x10.Hugo_Symbol,'verticalalignment','top','fontsize',9)

for i=1:nc

    x0=x10.cancer_cell_frac(i);
    y0=0*x10.cancer_cell_frac(i);

    x1=cosd(a);
    x2=0*x1;
    k=find((x1>0));
    x2(k)=x1(k).*(x10.high_CCF_CI(i)-x0);
    k=find((x1<=0));
    x2(k)=x1(k).*(x0-x10.low_CCF_CI(i));

    y1=sind(a);
    y2=y1*0.01;
    patch(x0+x2,y0+y2,cm(i,:),'edgecolor','none');
end

nc=slength(x20);
cm=jet(nc);

for i=1:size(x20.cancer_cell_frac)
vals3{i}=sprintf('%d-0',x20.cancer_cell_frac(i));
end
[cs,strs]=count(vals3);

for i=1:size(strs)
   color_vector3(ismember(vals3,strs{i}))=cs(i);
end

Xscat=vertcat(x12.cancer_cell_frac,x10.cancer_cell_frac,0*x20.cancer_cell_frac);
Yscat=vertcat(x21.cancer_cell_frac,0*x10.cancer_cell_frac,x20.cancer_cell_frac);
Cscat=vertcat(color_vector1',color_vector2',color_vector3');
Cs=(Cscat>1);

%plot(0*x20.cancer_cell_frac,x20.cancer_cell_frac,'k','.')
%text(0*x20.cancer_cell_frac+0.02,x20.cancer_cell_frac+0.02,x20.Hugo_Symbol,'verticalalignment','top','fontsize',9)
for i=1:nc

    x0=0*x20.cancer_cell_frac(i);
    y0=x20.cancer_cell_frac(i);

    x1=cosd(a);
    x2=0*x1;
    x2=x1*0.01;

    y1=sind(a);
    y2=0*y1;
    k=find((y1>0));
    y2(k)=y1(k).*(x20.high_CCF_CI(i)-y0);
    k=find((y1<=0));
    y2(k)=y1(k).*(y0-x20.low_CCF_CI(i));
    patch(x0+x2,y0+y2,cm(i,:),'edgecolor','none');
end
scatter(Xscat,Yscat,60,Cs,'fill')
colormap([0.75,0.75,0.75; 0,0,0])
hPatch = findobj(gca,'Type','patch');
set(hPatch,'facealpha',0.08)
axis equal;
axis([-.01 1.1 -.01 1])
title(pair_mapping_file.individual_id{zz},'FontSize',20)
xlabel(pair_mapping_file.pair_id_p{zz},'FontSize',20)
ylabel(pair_mapping_file.pair_id_m{zz},'FontSize',20)
f=sprintf('/xchip/cga_home/amaro/EwingsSarcoma/%s_EwingsCloud.png',pair_mapping_file.individual_id{zz});
%print_D(f,{{'fig'},{'pdf'},{'png','-r300'}});

%saveas(gcf,f,'eps2c')
print(gcf, '-dpng', '-r400',[f,'.png'])
% frame=getframe(gcf);
% if isempty(frame.colormap)
%    imwrite(frame.cdata, f)
% else
%    imwrite(frame.cdata, frame.colormap, f,'Quality',100)
% end

hold off
% 
% end
end
% 
% %{%% MY OWN PLOT:
% z=repmat('|',length(X0.Chromosome),1);
% X0.site=regexprep(cellstr([char(X0.Chromosome) z num2str(X0.Start_position) z char(X0.Tumor_Seq_Allele1)]),' ','');
% pri_ind= ismember(X0.patient, SETD.pair_id);
% met_ind = ismember(X0.patient, SETT.pair_id);
% pri = reorder_struct(X0, pri_ind);
% met = reorder_struct(X0, met_ind);
% 
% [o met_ind xom_ind] = intersect(met.site, X0.site);
% [o pri_ind xop_ind] = intersect(pri.site, X0.site);
% 
% total=struct;
% total.site = X0.site;
% total.Hugo_Symbol = X0.Hugo_Symbol
% total.CCF_met = zeros(length(total.site),1);
% total.CCF_pri = zeros(length(total.site),1);
% total.CCF_met(xom_ind) = met.CCF(met_ind);
% total.CCF_pri(xop_ind) = pri.CCF(pri_ind);
% 
% %Heirarchical clustering:
% X = [total.CCF_pri total.CCF_met];
% %Z = linkage(X,'single');
% %T = cluster(Z, 'maxclust', 4);
% 
% seed = [0 0; 0 1; 1 0; 1 1];
% [T cntr] = kmeans(X,4,'distance','cityblock','start',seed);
% X0.clusters = T;
% total.clusters = T;
% 
% clusts = unique(T);
% col = {'ko','mo','bo','ro'}
% 
% for i=1:length(clusts)
%     p = reorder_struct(total, total.clusters==clusts(i));
%     plot(p.CCF_pri, p.CCF_met, col{i});
%     hold on
% end
% 
% xlabel('Primary CCF')
% ylabel('Met CCF')
% 


%}

%X0=mergeStruct(XD,XT);

%X0.multiplicity=NaN*zeros(size(X0.id));
%X0.subclonal_ix=NaN*zeros(size(X0.id));
%X0.q_hat=NaN*zeros(size(X0.id));
%X0.modal_q_s=NaN*zeros(size(X0.id));

%[i m]=ismember(X0.id,X.id);
%k=find(i);
%X0.subclonal_ix(k)=X.subclonal_ix(m(m>0));
%X0.purity(k)=X.purity(m(m>0));
%X0.q_hat(k)=X.q_hat(m(m>0));
%X0.modal_q_s(k)=X.modal_q_s(m(m>0));
%X0.cancer_cell_frac(k)=X.cancer_cell_frac(m(m>0));
%X0.homozygous_ix(k)=X.homozygous_ix(m(m>0));
%X0.multiplicity=X0.cancer_cell_frac.*X0.modal_q_s;

 %X0.AF_Tumor=X0.i_tumor_f./X0.purity;
  
%X0.multiplicity(k)=X0.i_tumor_f(k).*X0.CN(k)./X0.purity(k);
%
%X0.homozygous_ix(k)=round(X0.multiplicity(k))==X0.CN(k) ;
%X0.modal_q_s(k)=round(X0.CN(k).*X0.purity(k))+ones(size(X0.CN(k)));
 %X0.modal_q_s(k)=round((X0.CN(k)-ones(size(X0.CN(k)))*2)./X0.purity(k) + 2*ones(size(X0.CN(k))));
 %X0.modal_q_s(X0.modal_q_s<1)=1;
%X0.cancer_cell_frac=X0.multiplicity./X0.modal_q_s(k);
%X0.cancer_cell_frac(find(X0.cancer_cell_frac>1))=1;
%and to convert to CCF we can use something like this with modal_q_s:
%X0.multiplicity=X0.cancer_cell_frac.*X0.modal_q_s;
