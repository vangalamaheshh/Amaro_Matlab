function s=trimatch(s,ch)

t=find(s==ch,1,'first');
if ~isempty(t)
    s=s(1:t-1);
end

