cd ~/Projects/TestesRotation/RSEM   
files = dir('.');
CDKN2A='ENSG00000147889.12';
POU5F1='ENSG00000204531.11';

for j=5:length(files)-3
  
r=load_struct(files(j).name);
i=j-4;
R.sample_id{i,1}=files(j).name;
R.POU5F1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000204531.11'))});
R.CDKN2A_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,CDKN2A))});              
R.PTEN_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000171862.5'))});
R.CDKN1A_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000124762.9'))});
R.TP53_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000141510.11'))});
R.VDR_RPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000111424.6'))});
R.NANOG(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000111704.6'))}); 
R.SOX2(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000181449.2'))}); 
R.SOX17(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000164736.5'))});
R.UTF1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000171794.3'))});
R.WNT1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000125084.7'))});
R.CTNNB1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000168036.12'))});


%Down
R.CCND2(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000118971.3'))});
R.WNT2(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000105989.4'))});
R.WNT2B(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000134245.13'))});

%up
Rup.TCF7L1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000152284.4'))});
Rup.FZD5_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000163251.3'))});
Rup.CXXC4_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000168772.9'))});
Rup.FZD8_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000177283.4'))});
Rup.FZD1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000157240.2'))});
Rup.TLE1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000196781.9'))});
Rup.FZD7_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000155760.1'))});
Rup.CTBP1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000159692.11'))});
Rup.LRP5_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000162337.7'))});
Rup.FZD4_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000174804.3'))});
Rup.SLC9A3R1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000109062.5'))});
Rup.SFRP1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000104332.7'))});
Rup.CTNNB1_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000168036.12'))});
Rup.WNT11_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000085741.8'))});
Rup.PITX2_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000164093.11'))});
Rup.TLE2_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000065717.10'))});
Rup.AES_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000104964.10'))});
Rup.PORCN_TPM(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000102312.16'))});








end
save_struct(R,'OCT34_CDKN2A_PTEN_TPM.txt')



%Up
R.TCF7L1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000152284.4'))});
R.FZD5(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000163251.3'))});
R.CXXC4(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000168772.9'))});
R.FZD8(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000177283.4'))});
R.FZD1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000157240.2'))});
R.TLE1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000196781.9'))});
R.FZD7(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000155760.1'))});
R.CTBP1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000159692.11'))});
R.LRP5(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000162337.7'))});
R.FZD4(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000174804.3'))});
R.SLC9A3R1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000109062.5'))});
R.SFRP1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000104332.7'))});
R.CTNNB1(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000168036.12'))});
R.WNT11(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000085741.8'))});
R.PITX2(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000164093.11'))});
R.TLE2(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000065717.10'))});
R.AES(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000104964.10'))});
R.PORCN(i,1)=str2num(r.TPM{find(ismember(r.gene_id,'ENSG00000102312.16'))});


















R=load_table('~/Downloads/OCT34_CDKN2A_PTEN_TPM.txt');
R=rmfield(R,{'headline','header'});
for i=1:slength(R)
strs=regexp(R.sample_id{i},'.rsem.genes.results','split');
id=strs{1};
R.sample_id{i}=id;
end

Stable=load_struct('~/Documents/GCT_Final_STable1.txt');

R_cohort=reorder_struct(R,ismember(R.sample_id,Stable.case_sample_id));

[i m]=ismember(R_cohort.sample_id,Stable.case_sample_id);

R_cohort.loc(i)=Stable.Testes_Sample(m);
R_cohort.histo(i)=Stable.Histology(m);
R_cohort.vital(i)=Stable.Vitalstatus(m);

boxplot(R_cohort.NANOG,R_cohort.loc)


boxplot(R_cohort.TP53_TPM,R_cohort.histo,'GroupOrder',{'Seminoma','Non-Seminoma mixed','Non-Seminoma pure','Teratoma','Other'})

RNA_MATRIX=load('/Volumes/xchip_cga_home/amaro/TestesRotation/RSEM/RNA_MATRIX.tpm.mat');
colids=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/RSEM/gene_col_ids');
rna=RNA_MATRIX.RNA_MATRIX;
testes_ind=find(ismember(R.sample_id,Stable.case_sample_id(ismember(Stable.Testes_Sample,'yes'))));
met_ind=find(ismember(R.sample_id,Stable.case_sample_id(ismember(Stable.Testes_Sample,'no'))));

for i=1:length(rna)
    summary_rna.gene{i,1}=colids.gene_id{i};
    summary_rna.tests(i,1)=mean(rna(testes_ind,i));
    summary_rna.mets(i,1)=mean(rna(met_ind,i));
    %summary_rna.ttest(i,1)=mwwtest(rna(testes_ind,i),rna(met_ind,i));

end


% boxplot(R_cohort.POU5F1_TPM,R_cohort.histo)
% 
% % 
% % ylabel('Transcripts per Million','FontSize',20)
% % xlabel('Histology','FontSize',20)
% boxplot(R_cohort.POU5F1_TPM,R_cohort.loc)
% ylabel('Transcripts per Million','FontSize',20)
% xlabel('Sample Location','FontSize',20)
% boxplot(R_cohort.CDKN2A_TPM,R_cohort.loc)
% ylabel('Transcripts per Million','FontSize',20)
% set(gca,'XTickLabel',{'Testes','Metastasis'},'FontSize',20)
% title('CDKN2A Expression','FontSize',30)
% boxplot(R_cohort.CDKN2A_TPM,R_cohort.histo)
% ylabel('Transcripts per Million','FontSize',20)
% xlabel('Histology','FontSize',20)
% title('CDKN2A Expression','FontSize',30)