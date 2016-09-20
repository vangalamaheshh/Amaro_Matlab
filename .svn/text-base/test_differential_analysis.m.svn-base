D.dat=rand(40000,20);

%tic;
%P1=differential_analysis(D,1:10,11:20,'ttest');
%t(1)=toc; % 73.3456
  
tic;
P2=differential_analysis_by_parts(D,1:10,11:20,'ttest',2,'/xchip/data/gadgetz/lsfres/');
t(2)=toc; 

tic;
P2=differential_analysis_by_parts(D,1:10,11:20,'ttest',4,'/xchip/data/gadgetz/lsfres/');
t(3)=toc; 

tic;
P2=differential_analysis_by_parts(D,1:10,11:20,'ttest',8,'/xchip/data/gadgetz/lsfres/');
t(4)=toc; 

tic;
P2=differential_analysis_by_parts(D,1:10,11:20,'ttest',16,'/xchip/data/gadgetz/lsfres/');
t(5)=toc;

tic;
P2=differential_analysis_by_parts(D,1:10,11:20,'ttest',32,'/xchip/data/gadgetz/lsfres/');
t(6)=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D.dat=rand(40000,20);
put_nan_at=find(rand(size(D.dat))>0.8);
D.dat(put_nan_at)=NaN;
t=[];
tic;
[PR,rs]=differential_analysis_permutations_by_parts(D,1:10,11:20,'ttest',500,1,'/xchip/data/gadgetz/lsfres/');
t(1)=toc


tic;
[PR,rs]=differential_analysis_permutations_by_parts(D,1:10,11:20,'ttest',500,10,'/xchip/data/gadgetz/lsfres/');
t(2)=toc


