function [Dobs P]=TajimasD(N,SNP_data,race)
%D calculation
S=slength(SNP_data);
pi=((2*N)/(N-1))*sum(SNP_data.(race).*(1.-SNP_data.(race)));
a1=sum(1./[1:N-1]);
a2=sum(1./([1:N-1]).^2);
b1=(N+1)/(3*(N-1));
b2=2*(N^2+N+3)/((9*N)*(N-1));
c1=b1-(1/a1);
c2=b2-((N+2)/(a1*N))+(a2/a1^2);
e1=c1/a1;
e2=c2/(a1^2+a2);
theta=S/a1;
Dobs=(pi-theta)/sqrt((e1*S)+e2*S*(S-1));
%P value calculation limit as S approaches inf
a=((2/N) - (1/a1))/e2^(1/2);
if mod(N,2)==0
   b= (N/(2*(N - 1)) - 1/a1)/e2^(1/2);
else
   b= ((N + 1)/(2*N) - 1/a1)/e2^(1/2);
end
alpha=-((1+a*b)*b)/(b-a);
Beta=((1+a*b)*a)/(b-a);
D=[a:.01:b];
BetaD=(gamma(alpha + Beta)*(b - D).^(alpha - 1).*(D - a).^(Beta - 1))/((gamma(alpha)*gamma(Beta)*(b - a)^(alpha + Beta - 1)));
cdfD=cumsum(BetaD./sum(BetaD));
if(Dobs < 0)
    P=cdfD(find(D<=Dobs,1,'last'));
else
   P=cdfD(find(D>=Dobs,1,'first'));
end
end
