function D1=discretize_D(D,params)

D1=D;
switch params.method
 case 'prctile'
  y=prctile(D.dat',params.prctile);
  D1.dat=nan(size(D.dat));
  D1.dat(D.dat<=repmat(y(1,:)',1,size(D.dat,2)))=0;
  D1.dat(D.dat>=repmat(y(2,:)',1,size(D.dat,2)))=1;
 case 'gap'
  y=prctile(D.dat',params.prctile);
  D1.dat=nan(size(D.dat));
  D1.dat(D.dat<=repmat(y(1,:)',1,size(D.dat,2)))=0;
  D1.dat(D.dat>=repmat(y(2,:)',1,size(D.dat,2)))=1;
  ltmat=repmat(y(1,:)',1,size(D.dat,2));
  gtmat=repmat(y(2,:)',1,size(D.dat,2));
  lt=find(D.dat<=ltmat);
  gt=find(D.dat>=gtmat);
  
  D.dat(lt)=ltmat(lt);
  D.dat(gt)=gtmat(gt);
  s=sort(D.dat,2);
  sn=(s-repmat(y(1,:)',1,size(s,2)))./repmat(diff(y,1)',1,size(s,2));
  
  for i=1:size(s,1)
    si=s(i,:);
    pos_hi=min(find(si==si(end)));
    pos_low=max(find(si==si(1)));
    xi=repmat(1:size(si,2),size(si,2),1);
    xj=xi';
    x=si(xi)-si(xj);
    d=xi-xj;
    sl1=x./(d+eps);
    repmat(si,size(si,2),1)-repmat(si',1,size(si,2));
%    imagesc(double(x>0 & x<params.gap*diff(y(:,i),1)));
    sl2=(si-y(1,i))./(pos_hi-(1:length(si))+eps);
    sl3=(si-y(2,i))./(pos_low-(1:length(si))+eps);
    
    figure(1); clf;
    subplot(2,2,1:2);
    imagesc(sl1);
    subplot(2,2,3);
    plot(1:(pos_hi-1),sl2(1:(pos_hi-1)));
    subplot(2,2,4);
    plot((pos_low+1):length(si),sl3((pos_low+1):end));
    
    
    n=sum(double(x<0 & x>-params.gap*diff(y(:,i),1)),2);
    n=sum(double(x>0 & x<params.gap*diff(y(:,i),1)),2);
    [m,mi]=min(n);
    th(i)=si(mi);
  end
 case 'slopes'
  keyboard
  sd=sort(D.dat,2);
  n=sum(~isnan(sd),2);
  E=zeros(size(D.dat,1),size(D.dat,2),size(D.dat,2));
  for s=1:size(D.dat,1)
    for i=(params.left+1):(n(s)-params.right-1)
      for j=(i+1):(n(s)-params.right)
        yleft=sd(s,1:(i-1));
        pleft=polyfit(1:(i-1),yleft,1);
        E(s,i,j)=E(s,i,j)+sum((yleft-pleft(1)*(1:(i-1))-pleft(2)).^2,2);
        ymid=sd(s,i:j); 
        pmid=polyfit(i:j,ymid,1);
        E(s,i,j)=E(s,i,j)+sum((ymid-pmid(1)*(i:j)-pmid(2)).^2,2);
        yright=sd(s,j:n(s)); 
        pright=polyfit(j:n(s),sd(s,j:n(s)),1);
        E(s,i,j)=E(s,i,j)+sum((yright-pright(1)*(j:n(s))- ...
                               pright(2)).^2,2);
      end
    end
    
    figure(1); clf;
    plot(1:n(s),sd(s,1:n(s)),'x'); hold on;
    plot(1:(i-1),pleft(1)*(1:(i-1))+pleft(2),'r-');
    plot(i:j,pmid(1)*(i:j)+pmid(2),'g-');
    plot(j:n(s),pright(1)*(j:n(s))+pright(2),'b-');
    disp([ s i j ]); pause;
  end
 case 'val'
  D1.dat=nan(size(D.dat));
  r=params.ranges;
  for i=1:size(r,1)
    D1.dat(find(D.dat>=r(i,1) & D.dat<=r(i,2)))=r(i,3);
  end
end
