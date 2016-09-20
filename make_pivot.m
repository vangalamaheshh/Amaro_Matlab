function [pivot,idx]=make_pivot(D,params)

switch params.method
 case 'by_type'
  if isfield(params,'supvals')
    v=params.supvals;
  else
    v=unique_keepord(D.supdat(params.supid,:));
    v=v(~isnan(v));
  end
  idx=find(ismember(D.supdat(params.supid,:),v));
  p=aggregate_D(reorder_D_cols(D,idx),'cols',params.supid,params.aggtype);
  if length(v)>1
    if isfield(params,'supaggtype')
      p=aggregate_D(p,'cols',ones(1,size(p.dat,2)),params.supaggtype);
    else
      p=aggregate_D(p,'cols',ones(1,size(p.dat,2)),params.aggtype);
    end
  end
  pivot=p; %p.dat;
 otherwise
  error('none');
end

