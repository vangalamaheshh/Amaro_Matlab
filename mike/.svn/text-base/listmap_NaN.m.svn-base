function m = listmap(a,b)

if ischar(a), a={a}; end
if ischar(b), b={b}; end

m = nan(length(b),1);
[a ai aj] = unique(a);
[c ia ib] = intersect(a,b);

for i=1:length(c), if ~mod(i,1e5), fprintf('%d/%d ',i,length(c)); end
  
    
    m(aj==ia(i)) = ib(i);
end, if length(c)>=1e5, fprintf('\n'); end


% for i=1:length(m)
%     if isnan(m(i))
%         m(i:end)=m(i:end)+1;
%     end
% end
m(isnan(m))=[];    




return

% old method
m = nan(length(a),1);
[u ui uj] = unique(a);
for i=1:length(u)
  idx = find(strcmp(u{i},b),1);
  if isempty(idx), idx = NaN; end
  m(uj==i) = idx;
end
