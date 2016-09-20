%% HWE function
function [p, chisq]=HWE_genetics(AA,AB,BB)

N=AA+AB+BB;
AFmajor=(AA+AA+AB)/(2*N);
eAA=AFmajor^2*(N);
eAB=2*AFmajor*(1-AFmajor)*N;
eBB=((1-AFmajor)^2)*N;
chisq=(eAA-AA)^2/eAA+(eAB-AB)^2/eAB+(eBB-BB)^2/eBB;
p=1-chi2cdf(chisq,1);

end
%% Code for Pset 5
function complete_problem_3

%% A
% determine HWE for SNPs
gc=load_table('~/Downloads/genotypes.counts.csv',','); % counts per genotype per snp
for i=1:slength(gc)
    gc.p(i,1)=HWE_genetics(gc.AA(i),gc.AB(i),gc.BB(i));
end
% divide by 100 to correct for the 100 multiple hypothesis so appropriate p value is
% .05/100
% gc.snp(find(gc.p<(.05*.01)))
% ans = 
%     'rs26'
%     'rs70'

%% B
gc2=load_table('~/Downloads/genotypes.csv',','); % genotypes per individual 
% Find individual with parents in genotype matrix
% parents should share ~50% of genotypes with individual 
% strategy: find rare SNPs where NA1021 has the minor allele 
gc.minor_af=(gc.BB+gc.AB+gc.BB)/50;
gcmaf=sort_struct(gc,'minor_af',1);
% gcmaf.snp(1:10)
% ans = 
%     'rs65'
%     'rs100'
%     'rs17'
%     'rs38'
%     'rs73'
%     'rs9'
%     'rs11'
%     'rs29'
%     'rs58'
%     'rs83'
% gc2.NA1021(ismember(gc2.snp,gcmaf.snp(1:10)))
% ans = 
%     'AA'
%     'CC'
%     'AG' <---- heterozygous rare alleles
%     'CC'
%     'CT' <---- heterozygous rare alleles
%     'AG' <---- heterozygous rare alleles
%     'AA'
%     'CC'
%     'CC'
%     'AA'
% find other patients with those rare alleles
SNPS=gc2.NA1021(ismember(gc2.snp,gcmaf.snp(1:20)));
idingSNPs=IDS(find(cellfun(@length,cellfun(@unique,SNPS,'UniformOutput',false))>1));
%NA1024 shares the rare allele @ rs17 also GG at rs40 also @ rs58
%NA1022 is the only other member who shares the allele  T at rs38 and also
%has the major allele at rs40. 
% snp	rs17	rs38	rs40	rs58
% NA1021	AG	CT	AG	AG
% NA1022	AA	CT	AA	AA
% NA1024	AG	CC	GG	AG

%% C
% find the SNPs that are in perfect LD (R^2 =1) with rs93
% correlation coefficient (SNPS,rs93)
rs_gc=load_table('~/Documents/rs_gc.txt');
fs=fieldnames(rs_gc);
fs=fs(2:end);

for i=1:length(fs)
    R=corrcoef(rs_gc.(fs{i}),rs_gc.rs93);
    CR(i,1)=R(1,2)^2;
end

% fs(CR==1)
% ans = 
%     'rs93'
%     'rs94'
%     'rs95'
%     'rs96'
%     'rs99'
% These snps are in perfect LD with rs93
end