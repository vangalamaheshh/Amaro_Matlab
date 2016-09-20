


%P=load_struct('Cranio_purity_file.txt')%('~/Projects/Cranios/Cranio_purity_file.txt');%Ewings PP call chart
P.purity=[0.52;0.28];
P.sample={'primary';'relapse'};
P.ploidy=[2;2];
P.Tumor_Sample_Barcode={'CNS-CN-001-Tumor-SM-7G9I4';'CNS-CN-001-Tumor-SM-7G9I5'};
P.maf={'~/Downloads/CNS-CN-001-TP-NB-SM-7G9I4-SM-7AABI.pon_filtered.txt';'~/Downloads/CNS-CN-001-TR-NB-SM-7G9I5-SM-7AABI.pon_filtered.txt'};
%SEGS=load_struct('~/Projects/Cranios/Exomes.segments.tsv');
%XD=load_table(''); %replace with Ewings Diagnostic MAF
%XT=load_table(''); %replace with relapse MAF
%X0=load_struct('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN/2830114/0.MutSigRun.Finished/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.final_analysis_set.maf');
p=load_struct(P.maf{1});
p2=load_struct(P.maf{2});
X0=load_struct('NEJ_Cranio_CCF_final');
X0=rmfields_if_exist(X0,'N');%('/xchip/cga/gdac-prod/cga/jobResults/MutSigRun/Exomes/4285440/Exomes.final_analysis_set.maf');
%X0.Start_position=cellfun(@str2num,X0.Start_position);
X0=rmfields_if_exist(X0,'headline');
X0=rmfields_if_exist(X0,'header');

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
X0.t_ref_count=str2double(X0.t_ref_count);

X0.t_alt_count=str2double(X0.t_alt_count);

CCF_range=0:0.01:1;
CCF_raw=0:0.01:100;

z=repmat('|',length(X0.Chromosome),1);
X0.id=regexprep(cellstr([char(X0.Chromosome) z char(X0.Start_position) z char(X0.Tumor_Sample_Barcode)]),' ','');

X0.purity=NaN*zeros(size(X0.id));

%X0.homozygous_ix=NaN*zeros(size(X0.id));
X0.cancer_cell_frac=NaN*zeros(size(X0.id));


k=find(isnan(X0.cancer_cell_frac));
X0.pat=X0.Tumor_Sample_Barcode;
for i=1:length(X0.pat)
X0.pat{i,1}=X0.pat{i}(1:3);
end
for i=1:length(P.sample)
P.pat{i,1}=P.sample{i}(1:3);
end
[i m]=ismember(X0.Tumor_Sample_Barcode,P.Tumor_Sample_Barcode);

%G0=trimStruct(gender_info_file,m(m>0));
%X0.gender=G0.gender;
X0.purity(ismember(X0.Tumor_Sample_Barcode,P.Tumor_Sample_Barcode{1}))=P.purity(1);
X0.purity(ismember(X0.Tumor_Sample_Barcode,P.Tumor_Sample_Barcode{2}))=P.purity(2);


% for i=1:length(X0.segvalue)
% X0.capseg_val(i,1)=(2^X0.segvalue(i)*2)/X0.purity(i);
% end



%X0.t_ref_count=cellfun(@str2num,X0.t_ref_count);
%X0.t_alt_count=cellfun(@str2num,X0.t_alt_count);
%X0.i_tumor_f=cellfun(@str2num,X0.i_tumor_f);


X0.O_tumor=X0.t_alt_count./(X0.t_alt_count+X0.t_ref_count);


%X0.CN=X0.capseg_val.*(X0.ploidy./2);

for i=1:length(X0.t_alt_count)
 if isequal(X0.Chromosome{i},'X')
CCFPdist1=betapdf((CCF_range).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
[val1,tmo1]=max(CCFPdist1);
           X0.CCF(i,1)=CCF_range(tmo1);
           [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
           n_dist=CCFPdist1./sum(CCFPdist1);
   X0.clone_dist(i,1)=sum(n_dist(80:100));
else
CCFPdist1=betapdf((CCF_range./2).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
[val1,tmo1]=max(CCFPdist1);
           X0.CCF(i,1)=CCF_range(tmo1);
           [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
           n_dist=CCFPdist1./sum(CCFPdist1);
   X0.clone_dist(i,1)=sum(n_dist(80:100));
end
end
%Xs=ismember(X0.Chromosome,'X');
%Males=ismember(X0.gender,'Male');
%MX=find((Xs+Males)==2);
%X0.CN(MX)=X0.CN(MX)./2;
% for i=1:length(X0.t_alt_count)
%     
%    if isequal(X0.Chromosome{i},'X') 
%            
%            CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%            [val1,tmo1]=max(CCFPdist1);
%            X0.CCF(i,1)=CCF_range(tmo1);
% [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%    else
%     
%    CCFPdist1=betapdf((CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%    [val1,tmo1]=max(CCFPdist1);
%    CCFPdist2=betapdf(2.*(CCF_range./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%    [val2,tmo2]=max(CCFPdist2);
%    if val1>val2 || X0.CN(i)<1.5
% [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%  X0.CCF(i,1)=CCF_range(tmo1);
%    else
% [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%  X0.CCF(i,1)=CCF_range(tmo1);
%    end
%    end
%    n_dist=CCFPdist1./sum(CCFPdist1);
%    X0.clone_dist(i,1)=sum(n_dist(80:100));
%    
%    
%    
%    
% end
% 
% 
% for i=1:length(X0.t_alt_count)
%     
%    if isequal(X0.Chromosome{i},'X') 
%            
%            CCFPdist1=betapdf((CCF_raw./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%            
%            CCFPdist1(101)=sum(CCFPdist1(101:end));
%             CCFPdist1(102:end)=[];
%               [val1,tmo1]=max(CCFPdist1);
% 
%    X0.CCF_raw(i,1)=CCF_raw(tmo1);
% [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%    else
%     
%    CCFPdist1=betapdf((CCF_raw./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%    CCFPdist1(101)=sum(CCFPdist1(101:end));
%    CCFPdist1(102:end)=[];
%    [val1,tmo1]=max(CCFPdist1);
%    CCFPdist2=betapdf(2.*(CCF_raw./X0.CN(i)).*X0.purity(i),X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
%    CCFPdist2(101)=sum(CCFPdist2(101:end));
%    CCFPdist2(102:end)=[];
%    [val2,tmo2]=max(CCFPdist2);
%    
%    if val1>val2 || X0.CN(i)<1.5
% [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%  X0.CCF_raw(i,1)=CCF_raw(tmo1);
%    else
% [X0.phat(i,1), X0.pci{i,1}] = binofit(X0.t_alt_count(i),(X0.t_alt_count(i)+X0.t_ref_count(i)),.05);
%  X0.CCF_raw(i,1)=CCF_raw(tmo1);
%    end
%    end
%    n_dist=CCFPdist1./sum(CCFPdist1);
%    X0.clone_dist_raw(i,1)=sum(n_dist(80:100));
%    
%    
% end
 
for i=1:length(X0.pci)
    X0.total_count(i,1)=X0.t_alt_count(i)+X0.t_ref_count(i);
    X0.af_CI95_low(i,1) = X0.t_alt_count(i)-(X0.pci{i}(1)*X0.total_count(i,1));
    X0.af_CI95_high(i,1) = X0.pci{i}(2)*X0.total_count(i,1)+X0.t_alt_count(i); 
    
    if X0.af_CI95_high(i)>X0.total_count(i)
        X0.af_CI95_high(i)=X0.total_count(i);
    end
    if isequal(X0.Chromosome{i},'X')
            CCFdist_low=betapdf((CCF_range).*X0.purity(i),X0.af_CI95_low(i)+1,X0.t_ref_count(i)+1);
    CCFdist_high=betapdf((CCF_range).*X0.purity(i),X0.af_CI95_high(i)+1,X0.t_ref_count(i)+1);
    [~,tmo_low]=max(CCFdist_low);
    [~,tmo_high]=max(CCFdist_high);
    X0.high_CCF_CI(i,1)=CCF_range(tmo_high);
    X0.low_CCF_CI(i,1)=CCF_range(tmo_low);
    else
    CCFdist_low=betapdf((CCF_range./2).*X0.purity(i),X0.af_CI95_low(i)+1,X0.t_ref_count(i)+1);
    CCFdist_high=betapdf((CCF_range./2).*X0.purity(i),X0.af_CI95_high(i)+1,X0.t_ref_count(i)+1);
    [~,tmo_low]=max(CCFdist_low);
    [~,tmo_high]=max(CCFdist_high);
    X0.high_CCF_CI(i,1)=CCF_range(tmo_high);
    X0.low_CCF_CI(i,1)=CCF_range(tmo_low);
    end
end




X0.subclonal_ix(k,1)=X0.clone_dist(k)<0.5;
save_struct(X0,'/Users/amaro/Documents/Exomes.final_analysis_mafCCF10_11.txt');%'/xchip/cga_home/amaro/Cranios/CCF_annotatedMaf.txt')


