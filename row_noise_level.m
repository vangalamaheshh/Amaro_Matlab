function [sig,Neff]=row_noise_level(x)


N=size(x,2);
p=[25 50 75];
T=10000;
X=random('norm',0,1,N,T);

Nhalf=floor(N/2);
for j=1:2
    if j==2
        X=X(1:Nhalf,:);
    end
    S=sort(X,1);
    D=diff(S,1,1);
    for i=1:length(p)
        y{i,j}=prctile(D,p(i));
        factor(i,j)=1./mean(y{i,j});
    end
end

for i=1:length(p)
    P(i,1)=(factor(i,1)-factor(i,2))/(N-Nhalf);
    P(i,2)=factor(i,1)-N*P(i,1);
end

X=x';

S=sort(X,1);
D=diff(S,1,1);
for i=1:length(p)
    Y{i}=prctile(D,p(i));
end

i1=2;
%i2=1;
%keyboard

%Neff=(P(i1,2)*Y{i1}-P(i2,2)*Y{i2})./(Y{i2}*P(i2,1)-Y{i1}*P(i1,1));
%Neff=Neff';

%sig=Y{i1}'*(P(i1,2)+P(i1,1)*N);

sig=Y{i1}'*factor(i1,1);

