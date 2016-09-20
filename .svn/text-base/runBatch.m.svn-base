cd /xchip/data/gadgetz/new_partial/packSS/main/

lsfdir='/xchip/data/gadgetz/lsfres/';

l=lsf(lsfdir);
nparts=20;
h=zeros(nparts,1);
for i=1:nparts
%  q=mod((i-1),2)+2;
%  iter=floor((i-1)/2)+1;
  q=2;
  iter=i;
  [l,h(i)]=bsub(l,{'t'},'threeWaves',{q,iter},'long');
end
[l,res]=wait(l); % wait for all  
cd /xchip/data/gadgetz/new_partial/packSS/resTmpSemi
unix('tar cvfz tt2.tar.gz tt2*.mat; rm tt2*.mat');

nparts=20;
h=zeros(nparts,1);
for i=1:nparts
%  q=mod((i-1),2)+2;
%  iter=floor((i-1)/2)+1;
  q=3;
  iter=i;
  [l,h(i)]=bsub(l,{'t'},'threeWaves',{q,iter},'long');
end
[l,res]=wait(l); % wait for all  
cd /xchip/data/gadgetz/new_partial/packSS/resTmpSemi
unix('tar cvfz tt3.tar.gz tt3*.mat; rm tt3*.mat');


