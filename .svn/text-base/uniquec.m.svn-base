function [u,ui,uj,uc]=uniquec(A,is_rows)

if nargin==2
  [u,ui,uj]=unique(A,is_rows);
else
  [u,ui,uj]=unique(A);
end

uc=histc(uj,1:max(uj));
[ucs,uci]=sort(-uc,1);

uc=uc(uci);
ui=ui(uci);

if exist('is_rows','var') && strcmp(is_rows,'rows')
    u=u(uci,:);
else
    u=u(uci);
end 
uj2=zeros(size(uj));
for i=1:length(uci)
  uj2(find(uj==uci(i)))=i;
end
uj=uj2;
