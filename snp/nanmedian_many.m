function y=nanmedian_many(x)
% calc median along dim=2

N=size(x,2);
n=sum(~isnan(x),2);
y=zeros(size(x,1),1);

is0=find(n==0);
if ~isempty(is0)
  y(is0)=NaN;
end

for i=1:size(x,2)
  isi=find(n==i);
  L=length(isi);
  if ~isempty(isi)
    if mod(i,2) % odd
      for j=1:ceil(i/2)
        [mv,mi]=min(x(isi,:),[],2);
        x(sub2ind(size(x),isi,mi))=NaN;
      end
      y(isi)=mv;
    else % even
      mv=zeros(length(isi),1);
      for j=1:(i/2+1)
        prev=mv;
        [mv,mi]=min(x(isi,:),[],2);
        x(sub2ind(size(x),isi,mi))=NaN;
      end
      y(isi)=(prev+mv)/2;
    end
  end
end

