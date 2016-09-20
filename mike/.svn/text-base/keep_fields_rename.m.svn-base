function S2 = keep_fields_rename(S,flds1,flds2)
% keep_fields_rename(S,flds1,flds2)
%
% given struct <S> and cell-array-of-strings <flds1>,
% returns struct <S2> which has only those fields specified, renamed to flds2
%
% Mike Lawrence 2010-09-14

if length(flds1)~=length(flds2), error('flds1 and flds2 should be same length'); end

S2=[];
for i=1:length(flds1)
  if isempty(S)
    f = [];
  else
    f = getfield(S,flds1{i});
  end
  S2=setfield(S2,flds2{i},f);
end
