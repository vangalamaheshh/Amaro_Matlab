function S1 = struct_assign(S1,idx,S2,flds)
% struct_assign(S1,idx,S2[,flds])
%
% S1(idx)=S2
%
% Mike Lawrence 2010

if slength(S2)~=length(idx), error('length(idx) ~= slength(S2)'); end
if max(idx)>slength(S1), error('max(idx) > slength(S1)'); end

f1 = fieldnames(S1);

if ~exist('flds','var')
  f2 = fieldnames(S2);
else
  f2 = flds;
end

for i=1:length(f2)
  fn = f2{i};
  s2 = getfield(S2,fn);
  if ismember(fn,f1)
    s1 = getfield(S1,fn);
  else  % create new field
    if iscell(s2), s1 = cell(slength(S1),1); else s1 = nan(slength(S1),1); end
  end
  s1(idx) = s2;
  S1 = setfield(S1,fn,s1);
end
