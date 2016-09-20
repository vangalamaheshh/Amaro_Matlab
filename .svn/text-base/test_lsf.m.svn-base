clear class
clear l h

set_verbose_level(10);
l=lsf('/xchip/data/gadgetz/lsfres/',[],1);

n=5;
h=zeros(n,1);
for i=1:n
%  h(i)=bsub(l,{'x1','x2'},'sort',{x(:,i)}); % [x1,x2]=sort(x(:,i));
  [l,h(i)]=bsub(l,{'v','d'},'eig(dna_dist(rand(1000,1000)))'); % [v,d]=eig(dna_dist(rand(1000,1000)));
end
h

[l,res]=wait(l); % wait for all


