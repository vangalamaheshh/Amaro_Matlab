function [thst_amp,thst_del]=run_permutations(C,t1s,t2s,nperm)

disp('using the smoothed data!');

l=lsf('/xchip/data/gadgetz/lsfres/');
nparts=10;
h=zeros(nparts,1);
for i=1:nparts
  [l,h(i)]=bsub(l,{'hst_amp','hst_del'},'glioma_perm',{C.smooth,floor(nperm/nparts),...
                      struct('method','log',...
                             't1s',t1s,...
                             't2s',t2s)});
end
[l,res]=wait(l); % wait for all

thst_amp=zeros(1001,length(t1s)); % 4 t1s and t2s
thst_del=zeros(1001,length(t2s));
for i=1:10
  thst_amp=thst_amp+res{i}.hst_amp;
  thst_del=thst_del+res{i}.hst_del;
end
