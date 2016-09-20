


P=load_struct('Cranio_purity_file.txt')%('~/Projects/Cranios/Cranio_purity_file.txt');%Ewings PP call chart
%SEGS=load_struct('~/Projects/Cranios/Exomes.segments.tsv');
%XD=load_table(''); %replace with Ewings Diagnostic MAF
%XT=load_table(''); %replace with relapse MAF
%X0=load_struct('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN/2830114/0.MutSigRun.Finished/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.final_analysis_set.maf');
X0=load_table('/Users/amaro/Downloads/Exomes.final_analysis_set.maf');%('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/Exomes/4285440/Exomes.final_analysis_set.maf');
%X0.Start_position=cellfun(@str2num,X0.Start_position);
X0=rmfield(X0,'headline');
X0=rmfield(X0,'header');

% for i=1:length(SEGS.individual_id)
% MUTS=ismember(X0.patient,SEGS.individual_id{i});
% MUTS=find(MUTS==1);
% SEG=load_struct(SEGS.capseg_segmentation_file{i});
% SEG.End=cellfun(@str2num,SEG.End);
% SEG.Start=cellfun(@str2num,SEG.Start);
% 
%     for j=1:length(MUTS)
%     pos=X0.Start_position(MUTS(j));
%     chr=X0.Chromosome{MUTS(j)};
%     locs=find(ismember(SEG.Chromosome,chr)==1);
%     ll=find(SEG.End(locs)>pos,1,'first');
%     X0.segvalue{MUTS(j),1}=SEG.Segment_Mean{ll};
%  
%     
%     end
% 
%     
% end

%X0.segvalue=cellfun(@str2num,X0.segvalue);


CCF_range=0:0.01:1;

z=repmat('|',length(X0.Chromosome),1);
X0.id=regexprep(cellstr([char(X0.Chromosome) z char(X0.Start_position) z char(X0.Tumor_Sample_Barcode)]),' ','');

X0.purity=NaN*zeros(size(X0.id));

%X0.homozygous_ix=NaN*zeros(size(X0.id));
X0.cancer_cell_frac=NaN*zeros(size(X0.id));


k=find(isnan(X0.cancer_cell_frac));
X0.pat=X0.patient;
for i=1:length(X0.pat)
X0.pat{i,1}=X0.pat{i}(1:3);
end
for i=1:length(P.sample)
P.pat{i,1}=P.sample{i}(1:3);
end
[i m]=ismember(X0.pat,P.pat);
P0=trimStruct(P,m(m>0));
%G0=trimStruct(gender_info_file,m(m>0));
%X0.gender=G0.gender;
P0.purity=cellfun(@str2num,P0.purity,'UniformOutput', false);
P0.purity=cell2mat(P0.purity);
X0.purity(k)=P0.purity(k);
P0.ploidy=cellfun(@str2num,P0.ploidy,'UniformOutput', false);
P0.ploidy=cell2mat(P0.ploidy);
X0.ploidy=P0.ploidy;
% for i=1:length(X0.segvalue)
% X0.capseg_val(i,1)=(2^X0.segvalue(i)*2)/X0.purity(i);
% end



%X0.t_ref_count=cellfun(@str2num,X0.t_ref_count);
%X0.t_alt_count=cellfun(@str2num,X0.t_alt_count);
%X0.i_tumor_f=cellfun(@str2num,X0.i_tumor_f);


X0.O_tumor=X0.t_alt_count./(X0.t_alt_count+X0.t_ref_count);


%X0.CN=X0.capseg_val.*(X0.ploidy./2);
X0.CN=X0.ploidy;
%Xs=ismember(X0.Chromosome,'X');
%Males=ismember(X0.gender,'Male');
%MX=find((Xs+Males)==2);
%X0.CN(MX)=X0.CN(MX)./2;
for i=1:length(X0.t_alt_count)
    
   if isequal(X0.Chromosome{i},'X') 
           
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
   n_dist=CCFPdist1./sum(CCFPdist1);
   X0.clone_dist(i,1)=sum(n_dist(80:100));
   
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




X0.subclonal_ix(k,1)=X0.clone_dist(k)<0.5;
save_struct(X0,'/Users/amaro/Documents/Exomes.final_analysis_mafCCF.txt');%'/xchip/cga_home/amaro/Cranios/CCF_annotatedMaf.txt')


