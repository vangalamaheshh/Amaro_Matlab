function [pv,stat]=ranksum_many(dat,cl1,cl2,tail);

has_nan=sum(isnan(dat(:,[cl1 cl2])),2);

pv=zeros(size(dat,1),1);
stat=zeros(size(dat,1),1);
for k=1:size(dat,1)
  dat1=dat(k,cl1);
  dat2=dat(k,cl2);
  if has_nan(k)
    dat1=dat1(~isnan(dat1));
    dat2=dat2(~isnan(dat2));
    if ~isempty(dat1) & ~isempty(dat2)
      pv(k)=my_ranksum(dat1,dat2,0.01,tail);
      stat(k)=median(dat1)-median(dat2);
    else
      pv(k)=NaN;
      stat(k)=NaN;
    end
  else
    pv(k)=my_ranksum(dat1,dat2,0.01,tail);
    stat(k)=median(dat1)-median(dat2);
  end
%  pv(k)=ranksum(dat1,da2,0.01);
end

