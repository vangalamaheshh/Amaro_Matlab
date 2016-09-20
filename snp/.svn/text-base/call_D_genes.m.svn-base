function [Z,Y]=call_D_genes(CL,genes,ts,cyto,rg,is_rg_collapsed)

if ~exist('is_rg_collapsed','var') || isempty(is_rg_collapsed)
  rgg=collapse_rg_to_genes(rg);
  rg=order_rg_by_pos(rgg);
end

[mx,mn,Y]=get_high_cutoffs(CL,cyto);

tvals{1}=[ repmat(ts(1),1,size(CL.dat,2)); mx+ts(1)];
tvals{2}=[ repmat(ts(2),1,size(CL.dat,2)); -(mn-ts(2))];

[Z,Y]=call_genes(CL,rg,genes,tvals);
