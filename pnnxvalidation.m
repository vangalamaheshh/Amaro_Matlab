function pnnxvalidation(gctfile,clsfile,outfile,fets,sigs,weights,use_median,dist_type)

% read gctfile
D=read_mit_gct_file(gctfile);

% read clsfile
D=read_mit_cls_file(D,clsfile);

N=size(D.dat,2);
n_cat=D_n_cat(D,1);

for c=1:D_n_cat(D,1)
  cls1=find(D.supdat(1,:)==c);
  cls0=find(D.supdat(1,:)~=c);
  xvmat=ones(N)-eye(N);
  
  % fet ranks
  
  % 
end

% fets
% sigs
% weights



