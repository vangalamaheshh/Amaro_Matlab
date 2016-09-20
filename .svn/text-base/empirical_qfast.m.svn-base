function q=empirical_qfast(p,b,d,drange)

q=p;
if exist('d','var')
  for i=1:length(p)
    nlp=-log(p(i));
    bin=find(drange(1:(end-1))<=nlp & drange(2:end)>nlp);
    if ~isempty(bin)
      q(i)=sum(d(bin:end));
    else
      error('no such bin');
    end
  end
else
  N=length(b);
  for i=1:length(p)
    q(i)=(length(find(b<=p(i)))+1)/(N+2);
  end
end
