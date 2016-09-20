function n=D_n_cls(D,supid)

pos=find(D.supacc(supid,:)=='/');
if isempty(pos)
  n=2
else
  n=length(pos)+1;
end

%u=unique(D.supdat(supid,:));
%n=sum(~isnan(u),2);

