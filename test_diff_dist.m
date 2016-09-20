% generate data from a gaussian N(0,1) with N samples
% sort the values
% get differences
% find min or some percentile (p) min: p=0
% what is the distribution of this value
% take most likely value


N=100;
p=50;
T=10000;
X=random('norm',0,1,N,T);

S=sort(X,1);
D=diff(S,1,1);

y=prctile(D,p);

hist(y,50);
mean(y)
std(y)

sig=y*38.4615;

sig2=row_noise_level(X');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% two gaussians

N=100;
p=50;
T=10000;

X=[random('norm',0,1,N*0.5,T); random('norm',10,1,N*0.5,T); ];

S=sort(X,1);
D=diff(S,1,1);

y=prctile(D,p);

hist(y,50);
mean(y)
std(y)

sig=y*38.4615;
sig2=row_noise_level(X');

