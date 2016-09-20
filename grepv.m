function [res,resi]=grepv(reg_exp,strs,res_is_idx)

if ischar(strs)
  strs=cellstr(strs);
end

resi=find(cellfun('isempty',cat(1,regexp(strs,reg_exp))));
res=strs(resi);

if nargout==0
  disp(catn(res));
end

if exist('res_is_idx','var') && res_is_idx
  [res,resi]=exchange_vars(res,resi);
end

