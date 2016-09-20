for i=1:length(samples)
PnR.sample{i,1}=samples{i};
PnR.nmisTiN(i,1)=sum(ismember(AML_maf_c.Tumor_Sample_Barcode,samples{i})&ismember(AML_maf_c.Variant_Classification,'Missense_Mutation'));
PnR.nonTiN(i,1)=sum(ismember(AML_maf_c.Tumor_Sample_Barcode,samples{i})&ismember(AML_maf_c.Variant_Classification,'Nonsense_Mutation'));
PnR.silentTiN(i,1)=sum(ismember(AML_maf_c.Tumor_Sample_Barcode,samples{i})&ismember(AML_maf_c.Variant_Classification,'Silent'));
PnR.spliceTiN(i,1)=sum(ismember(AML_maf_c.Tumor_Sample_Barcode,samples{i})&ismember(AML_maf_c.Variant_Classification,'Splice_Site'));



PnR.nmisWASHU(i,1)=sum(ismember(WASHU.Tumor_Sample_Barcode,samples{i})&ismember(WASHU.Variant_Classification,'Missense_Mutation'));
PnR.nonWASHU(i,1)=sum(ismember(WASHU.Tumor_Sample_Barcode,samples{i})&ismember(WASHU.Variant_Classification,'Nonsense_Mutation'));
PnR.silentWASHU(i,1)=sum(ismember(WASHU.Tumor_Sample_Barcode,samples{i})&ismember(WASHU.Variant_Classification,'Silent'));
PnR.spliceWASHU(i,1)=sum(ismember(WASHU.Tumor_Sample_Barcode,samples{i})&ismember(WASHU.Variant_Classification,'Splice_Site'));

end