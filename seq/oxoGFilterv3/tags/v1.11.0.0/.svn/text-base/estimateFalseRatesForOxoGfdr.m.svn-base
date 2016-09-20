function [FP,FPci,FN,FNci] =  estimateFalseRatesForOxoGfdr(NALT,iART,fdrStruct,Noxo,NoxoCI,PoxoG)
% TODO: Update documentation
% function [] =  OxoGfdr(NALT,NART,iART,FDR1,PoxoG,Noxo,NoxoCI)
%   OxoGfdr removes (cut) mutations according to the 
%   'fdr' OxoG filter method:
%   Given a list of NALT (ALT count reads), NART (ALT count in OxoG orientaiton),
%   iART flags (C>A or G>T SNVs) OxoGfdr calculates the p-values (pox) that a given event is 
%   consistent with OxoG artifact according to PoxoG, the binomial expected ratio for NART/NALT.
%   The Benjamini-Hochberg procedure estimates the multiple test corrected false detection 
%   rate (qox) for a given pox threshold. The filter removes the highest N events with qox>FDR1
%
%   [Ncut,cut,pox,qox] =  OxoGfdr(NALT,NART,iART,FDR1,PoxoG) 
%   inputs: 
%       NALT: vector of ALT count reads in for each event
%       NART: vector ALT count reads in OxoG orientaiton for each event
%       iART: flags C>A or G>T SNVs
%       FDR1: target false detection rate for OxoG 'detection'. defaults to 0.01
%       PoxoG: expected ratio for NART/NALT
%       Noxo: estimated # OxoG from binomial mixture
%       NoxoCI: 95% CI # OxoG from binomial mixture
%
%   outputs:
%       cut: vector of flags that a mutation is cut (filtered out)
%       Ncut: number of cut events (optional)
%       pox: p-values that a given event is consistent with OxoG artifact
%       qox: Benjamini-Hochberg corrected for a given pox threshold 
%       FP:  Rejected Real calls
%       FPci:  95% Confidence interval for Rejected Real calls
%       FN:  Passed OxoG calls
%       FNci:  95% Confidence interval Passed OxoG calls
%

% if (nargin<4)
%     FDR1=0.01;
% end
% if (nargin<5)
%     PoxoG=0.96;
% end
% 
% pox=binocdf(NART,NALT,PoxoG);
% %pox(~iART)=0; 
% qox=calc_fdr_value(pox.*iART);
% cut=qox>FDR1;
% Ncut=sum(cut);
% 
% Ntot=length(NALT);
% Nmut=Ntot-Ncut;

Ncut = sum(fdrStruct.cut);
pox = fdrStruct.pox;
cut = fdrStruct.cut;
%p cutline
if all(cut == 0)
    pmin = 1;
else
    pmin=min(pox(cut));
end
if (Ncut<1), pmin=1; end


% number of real SNVs in artifact mode
NR=sum(iART)-Noxo;
acb=0:100;

% allele count distribution in the non-artifact mode
nac=hist(NALT(~iART),acb);

% Initialize NcutN, proportion of non-artifact mode mutations that are
%   under the cut line.  I.e. would be cut if in the artifact mode.
NcutN=0;
for a=acb
    if nac(a+1)<1,continue,end
    fac=binopdf(0:a,a,0.5);
    pac=binocdf(0:a,a,PoxoG);
    NcutN=NcutN+nac(a+1)*sum(fac.*(pac>pmin));
end

FP=NR*(NcutN/sum(nac));

% Do not allow a negative FP.
FP = max([FP 0]);

% OxoG calls passed
FN=Noxo-Ncut+FP;
% floor at zero (can't have <0 estimated FN)
FN=max([0 FN]);

% Since error is modeled as a poisson, use this to calculate the FP
%   confidence interval.
[~,FPci] = poissfit(FP);
if (FPci(2)>Ncut),FPci(2)=Ncut;end
if FP > Ncut, FP = Ncut; end

dLow=sqrt( (Noxo-NoxoCI(1))^2 + (FP-FPci(1))^2 );
dHigh=sqrt( (NoxoCI(2)-Noxo)^2 + (FPci(1)-FP)^2 );
FNci=[FN-dLow FN+dHigh];
if (FNci(1)<0),FNci(1)=0;end 
NT=length(NALT)-Noxo;
if (FNci(2)>NT),FNci(2)=NT;end 
    