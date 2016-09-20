function st=remove_quotes(st,firstlast)
if nargin==1
  firstlast=0;
end

if iscell(st)
  for i=1:length(st)
    st{i}=remove_quotes_one(st{i},firstlast);
  end
else
  st=remove_quotes_one(st,firstlast);
end


function st=remove_quotes_one(st,firstlast)

if firstlast && length(st>=2) && st(1)=='"' && st(end)=='"'
  st=st(2:(end-1));
else
  pos=find(st=='"');
  if ~isempty(pos)
    st=st(setdiff(1:length(st),pos));
  end
end

