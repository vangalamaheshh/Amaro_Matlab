function [idx,ph]=plot_top_corr(D,supid,n)

if length(supid)>1
  v=supid;
else
  v=D.supdat(supid,:);
end

cnv=dna_norm(v,1);
c=abs(dna_norm(D.dat,1)*cnv');

[cs,idx]=sort(1-c);

ph=plot(1:n,1-cs(1:n));
idx=idx(1:n);
