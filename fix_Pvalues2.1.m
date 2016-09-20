function [P1sgte,P1slte,P2s]=fix_Pvalues2(S,SR,smooth_flag,trueperm_is_first);

if (nargin>3) && ~trueperm_is_first
  SR=[S SR]; % true perm is first now
end

%if nnz(isnan(SR))
%  error('There was an NaN in the SR matrix');
%end

P1sgte=NaN*ones(size(SR));
P1slte=P1sgte;
P2s=P1sgte;

for i=1:size(P1sgte,1)
  not_nan=find(~isnan(SR(i,:)));
  Ns=length(not_nan);
  if not_nan>=1
    [u,a,b]=unique(SR(i,not_nan)'); % u is sorted and unique
    [pr,ord]=sort(b);
    [dum,revord]=sort(ord);
    n=length(u);
    h=diff(find(diff([0 pr' n+1])));
    pv=cumsum([0 h]); % ends with 1
    ipv=Ns-pv;
    pv=pv/Ns;
    ipv=ipv/Ns;
    tmp=ipv(1:(end-1));
    tmp2=other_tail(pv,ipv);
    P1sgte(i,not_nan)=tmp(b);
    tmp=pv(2:end);
    P1slte(i,not_nan)=tmp(b);
    P2s(i,not_nan)=min([P1sgte(i,not_nan); P1slte(i,not_nan)])+ ...
        tmp2(b)';
  end
end

