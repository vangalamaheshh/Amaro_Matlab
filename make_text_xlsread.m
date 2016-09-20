function r=make_text_xlsread(r)

emptypos=find(cellfun('isempty',r));
if ~isempty(emptypos)
  r(emptypos)=mat2cell(repmat(NaN,length(emptypos),1),ones(length(emptypos),1),1);
end

pos=find(cell2mat(cellfun_any('isnumeric(x)',r)));

n=cat(1,r{pos});
posnan=find(isnan(n));
r(pos(posnan))=cell(length(posnan),1);
pos_notnan=find(~isnan(n));
r(pos(pos_notnan))=cellstr(num2str(n(pos_notnan)));
