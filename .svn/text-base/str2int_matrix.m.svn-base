function n=str2int_matrix(x)

mask=ismember(double(x),(0:9)+48);

nonum=find(sum(mask,2)==0);
if ~isempty(nonum)
  warning([ num2str(length(nonum)) ' row(s) have no numbers']);
end

isnum=setdiff(1:size(x,1),nonum);
orig_x_sz=size(x,1);
x=x(isnum,:);
mask=mask(isnum,:);

x=double(x)-48;
tmp=diff([ zeros(size(x,1),1) mask zeros(size(x,1),1)],1,2);
[fi,fj]=ind2sub(size(tmp),find(tmp==1));
[sfi,sfi_idx]=sort(fi);
fj=fj(sfi_idx);
fi=sfi;

[li,lj]=ind2sub(size(tmp),find(tmp(:,2:end)==-1));
[sli,sli_idx]=sort(li);
lj=lj(sli_idx);
li=sli;

szi=lj-fj+1;
n=zeros(size(x,1),1);
cur=1:size(x,1);

p=1;
pos=sub2ind(size(x),li,lj);

while ~isempty(cur)
  n(cur)=n(cur)+p*x(pos(cur));
  p=p*10;
  szi=szi-1;
  cur=find(szi>0);
  pos=pos-size(x,1);
end


if ~isempty(nonum)
  nfull=nan(orig_x_sz,1);
  nfull(isnum)=n;
  n=nfull;
end
