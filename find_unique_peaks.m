function p=find_unique_peaks(x,wsz,cutoff)

p=zeros(size(x));

maxp=length(x)-wsz+1;
[mx,mxi]=max(x);
while mx>cutoff
  p(mxi)=1;
  
  if mxi<wsz || mxi>maxp
    x(clip_to_range((mxi-wsz+1):(mxi+wsz-1),[1 length(x)]))=0;
  else
    x((mxi-wsz+1):(mxi+wsz-1))=0;
  end
  [mx,mxi]=max(x);
end
