function [PoxoG,PoxoGci]=PoxoGestimate(NALT,NART,iART,alphaci)
%function [PoxoG,PoxoGci]=PoxoGestimate(NALT,NALT,iART)
%
% estimates underlying FoxoG probabilty (PoxoG)
%
% inputs:   NALT number of alternate alleles
%           NART number of alternate alleles in OxoG configuration
%           iART vector same size of NALT and NART =true for artifact-mode
%               (C>A, G>T), default everything is artifact mode
%           alphaci confidence interval alpha level (default 0.05)
% outputs:  PoxoG underlying OxoG bias toward F1R2 for G>T, F2R1 for C>A
%           PoxoGci confidence interval for PoxoG
%
if (nargin<3)
    iART=true(size(NART));
end
if (nargin<4)
    alphaci=0.05;
end

if sum(iART) == 0
   disp('No artifact mode mutations.  Cannot estimate PoxoG') 
   PoxoG = NaN;
   PoxoGci = NaN;
   return
end

[phat,pci] = binofit(NART(iART),NALT(iART),alphaci);
eART=diff(pci,[],2);
xART=NART(iART)./NALT(iART);
NL=sum(xART<0.5)+0.5*sum(xART==0.5); NL=max([NL 1]);
N=length(xART);
top=sum(xART)-NL;
sigma_top=sqrt(sum(eART.^2)+NL);
bot=N-2*NL;
sigma_bot=sqrt(N+2*NL);
PoxoG=top/bot;
PoxoGci = sqrt( (sigma_top/bot).^2 + (top*sigma_bot/bot^2).^2);

function test
f1='/local/cga-fh/cga/Test_OxoG/Individual_Set/Test13Jul2012/Individual/CESC-HSCX1127/jobs/capture/mut/oxog/filtered_maf/CESC-HSCX1127.oxoG.filter.maf.annotated'
f2='/local/cga-fh/cga/Test_OxoG/Individual_Set/Test13Jul2012/Individual/CESC-HSCX1127/jobs/capture/mut/oxog/filtered_maf/CESC-HSCX1127.oxoG.filter.maf.annotated.reject.maf.annotated'
X1=load_table(f1);
X2=load_table(f2);
X=mergeStruct(X1,X2)
X=trimStruct(X,find(ismember(X.Variant_Type,'SNP')))

TOFROM=strcat(X.Reference_Allele,X.Tumor_Seq_Allele1);
iCA=ismember(TOFROM,{'CA'});
iGT=ismember(TOFROM,{'GT'});
iART=ismember(TOFROM,{'CA','GT'});
NART=X.i_t_ALT_F1R2.*iGT + X.i_t_ALT_F2R1.*iCA;
NALT=X.i_t_ALT_F1R2+ X.i_t_ALT_F2R1;
alphaci=0.1;
[PoxoG,PoxoGci]=PoxoGestimate(NALT,NART,iART,alphaci)


