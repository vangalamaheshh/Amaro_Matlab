function rcmd=R_fisher_exact(x,y)

if ~(size(x,1)>1 && size(x,2)>1) % x is not a mat
  x=crosstab_full(x,y,1:max(max(x),max(y)));
end

s=max(size(x));
xf=zeros(s);
xf(1:size(x,1),1:size(x,2))=x;

rcmd=[ 'fisher.test(' R_mat(xf) ',workspace=10e6)' ];

