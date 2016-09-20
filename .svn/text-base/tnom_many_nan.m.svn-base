function [p,nerr]=tnom_many_nan(dat,cls0,cls1);

[sdat,idx]=sort(dat,2);
n0=sum(~isnan(dat(:,cls0)),2);
n1=sum(~isnan(dat(:,cls1)),2);

v=zeros(1,size(dat,2));
v(cls1)=1;
V=repmat(v,size(dat,1),1);

V=v(idx);

cs=[zeros(size(V,1),1) cumsum(V,2)];
err_np=2*cs+repmat(n0,1,size(dat,2)+1)-repmat(0:size(dat,2),size(cs,1),1);
err_pn=repmat(n1,1,size(dat,2)+1)+repmat(0:size(dat,2),size(cs,1),1)-2*cs;
err(:,:,1)=err_np;
err(:,:,2)=err_pn;
nerr=min(min(err,[],3),[],2);


if length(unique(n0))==1 && length(unique(n1))==1
  [u,a,b]=unique(nerr);
  for i=1:length(u)
    p(i)=calc_tnom_pv_closed2(n0(1),n1(1),u(i));
  end
  p=p(b)';
else
  p=zeros(size(dat,1),1);
  for i=1:size(p,1)
    p(i)=calc_tnom_pv_closed2(n0(i),n1(i),nerr(i));
  end
end

