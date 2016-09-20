function stgc=gc_std(mu,std0)
% took from differential_analysis.m and there from JPh
minstgc=0.2*abs(mu);
minstgc(minstgc<std0)=std0(minstgc<std0);
minstgc(minstgc==0)=0.1;

stgc=max([std0 minstgc],[],2);
