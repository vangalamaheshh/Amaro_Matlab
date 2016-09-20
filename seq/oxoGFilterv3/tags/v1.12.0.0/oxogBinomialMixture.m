function [Noxo,Nmut,alf,alfci,NoxoCI, lod0] =  oxogBinomialMixture(NALT,NART,AC,PoxoG,p0,alphaci)
% function [Noxo,Nmut,alf,alfci,eNoxo,lod0] =  oxogBinomialMixture(NALT,NART,AC,PoxoG,p0)
%   OxogBinomialMixture estimates the mixture of OxoG artifact and non-artifact 
%   components in a list of mutations with NALT ALT allele count reads and NART 
%   ALT OxoG orientation reads. This is calculated at ALT allele count
%   values listed in AC (default 1:100), for an OxoG expected NART/NALT
%   ratio of PoxoG (default 0.96), and a non-oxoG ratio of p0 (default 0.5)
%
%   [Noxo,Nmut,alf,alfci,eNoxo] =  OxogBinomialMixture(NALT,NART,AC,PoxoG,p0)
%   inputs: 
%       NALT: vector of ALT count reads for each event
%       NART: vector ALT count reads in OxoG orientaiton for each event
%       AC: vector of allele counts (default 1:100)
%       PoxoG: OxoG expected ratio for NART/NALT (default 0.96)
%       p0: non-OxoG expected ratio for NART/NALT (default 0.50)
%       alphaci: significance level for confidence intervals (default 0.05 for 95% CI)
%
%   outputs:
%       Noxo: fit number of OxoG SNVs in NALT,NART list
%       Nmut: fit number of non-OxoG SNVs: Nmut = length(NALT)-Noxo
%       alf: fit ratio of Oxo component for each value of AC
%       alfci: fit 1-sigma confidence ratio of Oxo component for each value of AC
%       NoxoCI: estimated uncertainty on Noxo from fit. 
%       lod0: log odds ratio corresponding to artifact at Noxo to no artifact at all 
%


if (nargin<3)
    AC=1:100;
end
if (nargin<4)
    PoxoG=0.96;
end
if (nargin<5)
    p0=0.5;
end
if (nargin<6)
    %alphaci=0.3173; %one sigma
    %alphaci=0.0455; %two sigma
    alphaci=0.05; % 95%
end

% drop in negative log likelihood for alphaci (natural log_e)
dLL=norminv(1-alphaci/2,0,1)^2/2;


alpha=0.0:0.001:1; 
Na=length(alpha);
NA=length(AC);

alf=NaN*zeros(NA,1);
alfci=NaN*zeros(NA,2);
Nm=zeros(NA,1);
cilim=log10(exp(1))/2;
Nox1=zeros(NA,1);
Nox1ci=zeros(NA,2);
vlow=zeros(NA,1);
vhigh=zeros(NA,1);

% initialize lod0
lod0=0;

for ac=1:NA
    k=find(NALT==ac);
    nox=0:ac; Nox=length(nox);
  
    if length(k)<1, continue; end
    NALT1=NALT(k);
    NART1=NART(k);
    N1=length(k);
    n1=binopdf(nox',ac,p0);
    n2=binopdf(nox',ac,PoxoG);

    % We would like to test to see if there 
    
    n=hist(NART1,0:ac,'log');
    NLL = @(a) sum(-log(poisspdf(n',N1*(1-a)*n1+N1*a*n2)));
    [alf(ac),fval0] = fminbnd(NLL,0,1);
    NLL0 = @(a) sum(-log(poisspdf(n',N1*(1-a)*n1+N1*a*n2)))-(fval0+dLL);
    
    % lod0(ac) is log like of fitted alpha - log like of alpha==0
    lod0=lod0-NLL(alf(ac))+NLL(0);
    
    if lod0 < 0
        disp(['lod0 is a negative number (' num2str(lod0) '). ac: ' num2str(ac) '   NLL(0): ' num2str(NLL(0)) '   -NLL(alf(ac)): ' num2str(-NLL(alf(ac)))  ]) 
    end
    
    clear fv
    av=0:0.1:1;
    for i=1:length(av)
        fv(i,1)=NLL0(av(i));
    end
    
    klim=find((fv)<100); klim=klim([1 end]);
    if (klim(1)>1)
        klim(1)=klim(1)-1;
    end
    if (klim(2)<length(fv))
        klim(2)=klim(2)+1;
    end
        
    if (fv(1)>0)
        [alfci(ac,1), fvalCI,ef] = fzero(NLL0,[av(klim(1)) alf(ac)]);
        if (ef~=1)
            fplot(NLL0,[0 1])
            keyboard
        end
        
    else
        alfci(ac,1)=0;
    end
    if (fv(end)>0)
        [alfci(ac,2), fvalCI,ef] = fzero(NLL0,[alf(ac) av(klim(2))]);
        if (ef~=1)
            fplot(NLL0,[0 1])
            keyboard
        end
    else
        alfci(ac,2)=1;
    end
    
    Nox1(ac)=alf(ac).*N1;  
    
    % estimate ci on N1
    e=alf(ac).*N1;
    e(e<1e-2)=1e-2;
    [q, N1ci]=poissfit(e,alphaci);
    dN1low=N1-N1ci(1);
    dN1high=N1ci(2)-N1;
    % alpha
    dalflow=alf(ac)-alfci(ac,1);
    dalfhigh=alfci(ac,2)-alf(ac);
    % propagate variance low and high separately (~ non-negative ci)
    vlow(ac,1)=N1^2*dalflow^2+alf(ac)^2*dN1low^2;
    vhigh(ac,1)=N1^2*dalfhigh^2+alf(ac)^2*dN1high^2;

    
end

Ntot=length(NALT);
Noxo=sum(Nox1);
Nmut=Ntot-Noxo;

eNoxoLow=sqrt(sum(vlow));
eNoxoHigh=sqrt(sum(vhigh));
% eNoxo=round(sqrt(sum(diff(Nox1ci,1,2).^2)));

NoxoCI=[Noxo-eNoxoLow Noxo+eNoxoHigh];
if (NoxoCI(1)<0), NoxoCI(1)=0;end;
if (NoxoCI(2)>Ntot), NoxoCI(2)=Ntot;end;


if lod0 < 0
    warning(['Final lod0 is a negative number.  Setting it to zero.'])
    lod0 = 0;
end
