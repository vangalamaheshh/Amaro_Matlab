M_full=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage_WhiteList.maf');
M=load_struct('/Volumes/xchip_cga_home/amaro/CLL/DeTiN_Figure/DeTiNMAF.maf');
sig_genes=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/1_09_2015_PoNCut4_WithSalvage/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt');
sig_genes=reorder_struct(sig_genes,str2double(sig_genes.q)<=.1);
sig_genes=reorder_struct(sig_genes,~ismember(sig_genes.gene,{'ADAM30';'NHS';'TTK'})); %genes which are not expressed / did not validate
samples_deTiN.sample=unique(M.Tumor_Sample_Barcode);
samples_MRD.sample=unique(M_full.Tumor_Sample_Barcode);
for i=1:slength(samples_MRD)
samples_MRD.is_GCLL{i,1}=regexp(samples_MRD.sample{i},'GCLL','match');
end
samples_MRD=reorder_struct(samples_MRD,~cellfun(@isempty,samples_MRD.is_GCLL));
for i=1:slength(samples_deTiN)
samples_deTiN.nMut(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i}));
samples_deTiN.nMut_mutect(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&ismember(M.i_failure_reasons,''));
samples_deTiN.nMut_deTiN(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&~ismember(M.i_failure_reasons,''));
samples_deTiN.ndbsnp(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&~ismember(M.dbSNP_RS,''));
samples_deTiN.notdbsnp(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&ismember(M.dbSNP_RS,''));


samples_deTiN.ndbsnp_pre(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&~ismember(M.dbSNP_RS,'')&ismember(M.i_failure_reasons,''));
samples_deTiN.ndbsnp_post(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&~ismember(M.dbSNP_RS,'')&~ismember(M.i_failure_reasons,''));
samples_deTiN.nDriver(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&ismember(M.Hugo_Symbol,sig_genes.gene)&ismember(M.Variant_Classification,coding)&ismember(M.Variant_Type,'SNP'));    
samples_deTiN.prenDriver(i,1)=sum(ismember(M.Tumor_Sample_Barcode,samples_deTiN.sample{i})&ismember(M.Hugo_Symbol,sig_genes.gene)&ismember(M.i_failure_reasons,'')&ismember(M.Variant_Classification,coding)&ismember(M.Variant_Type,'SNP'));    

end


for i=1:slength(samples_MRD)
    
samples_MRD.ndbsnp(i,1)=sum(ismember(M_full.Tumor_Sample_Barcode,samples_MRD.sample{i})&~ismember(M_full.dbSNP_RS,''));
samples_MRD.nMut(i,1)=sum(ismember(M_full.Tumor_Sample_Barcode,samples_MRD.sample{i}));  
samples_MRD.nDriver(i,1)=sum(ismember(M_full.Tumor_Sample_Barcode,samples_MRD.sample{i})&ismember(M_full.Hugo_Symbol,sig_genes.gene)&ismember(M_full.Variant_Classification,coding)&ismember(M_full.Variant_Type,'SNP'));    
    
end

distribution_ordered_plot(samples_deTiN.ndbsnp_pre./samples_deTiN.nMut_mutect,(samples_deTiN.ndbsnp_pre+samples_deTiN.ndbsnp_post)./samples_deTiN.nMut,samples_MRD.ndbsnp./samples_MRD.nMut)

set(gca,'XTick',[1,2,3],'XTickLabel',{'pre-DeTiN','post-DeTiN','MRD-Samples'},'YTick',[0,.25,.5])
ylabel('Rate of dbSNP per mut')

boxplot([samples_deTiN.prenDriver;samples_deTiN.nDriver;samples_MRD.nDriver],[ones(slength(samples_deTiN),1);2*ones(length(samples_deTiN.nDriver),1);3*ones(slength(samples_MRD),1)])
