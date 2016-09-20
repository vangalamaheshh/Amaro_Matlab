function make_mutation_spectrum_plot(Maf,order)


X1=load_struct(Maf);
code=get_coding_class_muts;
coding_X1=reorder_struct(X1,ismember(X1.Variant_Classification,code));
coding_X1_SNP=reorder_struct(coding_X1,ismember(coding_X1.Variant_Type,'SNP'));
coding_X1_SNP.mut=strcat(coding_X1_SNP.Reference_Allele,'->',coding_X1_SNP.Tumor_Seq_Allele2);
coding_X1_SNP.mutix(ismember(coding_X1_SNP.mut,'A->C')|ismember(coding_X1_SNP.mut,'T->G'))=1;
coding_X1_SNP.mutix(ismember(coding_X1_SNP.mut,'A->G')|ismember(coding_X1_SNP.mut,'T->C'))=2;
coding_X1_SNP.mutix(ismember(coding_X1_SNP.mut,'A->T')|ismember(coding_X1_SNP.mut,'T->A'))=3;
coding_X1_SNP.mutix(ismember(coding_X1_SNP.mut,'G->C')|ismember(coding_X1_SNP.mut,'C->G'))=4;
coding_X1_SNP.mutix(ismember(coding_X1_SNP.mut,'G->A')|ismember(coding_X1_SNP.mut,'C->T'))=5;
coding_X1_SNP.mutix(ismember(coding_X1_SNP.mut,'G->T')|ismember(coding_X1_SNP.mut,'C->A'))=6;
%6=Turquoise
%5=Yellow
%4=Red
%3=Purple
%2=Green
%1=Blue

%order=load_struct('~/Documents/CLL8_sample_order.tsv');
order=load_struct(order);
X1=reorder_struct(X1,ismember(X1.sample,order.sample));

samples_mut_ix=order.sample;
for i=1:length(samples_mut_ix)
for j=1:6
mut_sums_matrix(i,j)=sum(ismember(coding_X1_SNP.sample,samples_mut_ix{i})&(coding_X1_SNP.mutix==j)')./sum(ismember(coding_X1_SNP.sample,samples_mut_ix{i}));
end


end
C=[0 0 1;0 1 0;153/255 0 153/255;1 0 0;1 1 0;0 1 1];
H=bar(mut_sums_matrix,'stacked')
for k=1:6
set(H(k),'facecolor',C(k,:),'edgecolor',C(k,:))
end

box off
set(gca,'XTick',[])
set(gca,'YTick',[])


total_coverage=28834635; % number from mutsig
%Mclonal_sublclonal=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/final_analysis_set.maf');
%X1.clonal_ix=str2double(X1.ccf_median_hack)>=.85;
%for i=1:length(samples_mut_ix)
%    clonal_subclonal_rates(i,1)=sum(ismember(X1.sample,samples_mut_ix{i})&(X1.clonal_ix==1))/28.834635;
%    clonal_subclonal_rates(i,2)=sum(ismember(X1.sample,samples_mut_ix{i})&(X1.clonal_ix==0))/28.834635;
%end

%bar(clonal_subclonal_rates,'stack')

end