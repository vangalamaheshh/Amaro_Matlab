function [mx,mn,Y]=get_high_cutoffs(C,cyto,params)

if ~exist('params','var') || isempty(params)
  params.n_thresh=500;
end

C=reorder_D_rows(C,find(C.chrn<=22));
C=add_cyto(C,cyto);
C.chrarmn=C.armn+C.chrn*2;
collapse_type='nanmedian';
Y=collapse_D(C,'chrarmn',collapse_type); % 'mode'
Y=reorder_D_rows(Y,find(Y.n_collapse>params.n_thresh));

mx=max(Y.dat,[],1);
mn=min(Y.dat,[],1);
