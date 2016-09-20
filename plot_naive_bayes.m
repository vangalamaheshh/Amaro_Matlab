function plot_naive_bayes(D,pred,tr,posterior,cls,LP_xi_giv_c,P_c,fig_type)


if exist('LP_xi_giv_c','var') && exist('P_c','var')
  if isempty(LP_xi_giv_c)
    if isempty(cls)
      for k=1:length(P_c)
        ck=find(D.supdat(1,:)==k);
        cls.mu{k}=mean(D.dat(:,ck),2); % +rand(size(D.dat,1),1)*3;
        cls.sig{k}=std(D.dat(:,ck),0,2);%+5;
      end
    end
    for k=1:length(P_c)
      for i=1:size(D.dat,1)
        for s=1:size(D.dat,2)
          LP_xi_giv_c(k,i,s)=0.5*(-log(2*pi*cls.sig{k}(i).^2)-((D.dat(i,s)-cls.mu{k}(i))/cls.sig{k}(i)).^2); 
          % LP_xi_giv_c(k,i,s)=log(normpdf(D.dat(i,s),cls.mu{k}(i),cls.sig{k}(i)));
        end
      end
    end
  end
  
  for k=1:length(P_c)
    LP_x_giv_c(k,:)=sum(squeeze(LP_xi_giv_c(k,:,:)),1);
  end
  P_c_giv_x=exp(LP_x_giv_c).*repmat(P_c,1,size(D.dat,2));
  P_x=sum(P_c_giv_x,1);  
  P_c_giv_x=P_c_giv_x./repmat(P_x,3,1);
  max(abs(range(posterior-P_c_giv_x)))
  posterior=P_c_giv_x;
  [dum,pred]=max(posterior,[],1);
  xx=LP_xi_giv_c-repmat(shiftdim(log(P_x)/size(D.dat,1),-1),[length(P_c) size(D.dat,1) 1])+...
     repmat(log(P_c)/size(D.dat,1),[1 size(D.dat,1) size(D.dat,2)]);
  xx=exp(xx);
  xmin=prctile(xx(:),5);
  xmax=prctile(xx(:),95);
%  xmin=min(xx(:));
%  xmax=max(xx(:));
end

nc=size(posterior,1);
switch fig_type
 case 'matrices'
  gr=make_subplotgrid(ones(1,nc),[0.1 0.1 2 0.1 ],ones(1,nc+1),[1 1 1 1 1],0.05,0.05);
  for k=1:nc
    subplotgrid(gr,3,k);
    imagesc(squeeze(xx(k,:,:)));
    caxis([xmin xmax]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    subplotgrid(gr,4,k);
    imagesc(P_c_giv_x(k,:)==max(P_c_giv_x,[],1));
    caxis([0 1]);
    whitered;
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    subplotgrid(gr,2,k);
    imagesc(P_c_giv_x(k,:));
    caxis([0 1]);
    whitered;
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    subplotgrid(gr,1,k);
    v=zeros(1,size(D.dat,2));
    v(find(D.supdat(1,:)==k))=1;
    imagesc(v);
    caxis([0 1]);
    whitered;
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
  end
  
 case 'lines'
  gr=make_subplotgrid([ 0.4 ones(1,nc) 0.1],[0.1 0.1 2/3 2/3 2/3 0.1 ],ones(1,nc+3),[1 1 1 1 1 1 1],0.05,0.05);
  D1=preprocess_D(D,'row_center_and_normalize');
  subplotgrid(gr,3:5,1);
  display_D_elem(gr,'gacc',D1,[],[],'vert',8);
%  subplotgrid(gr,3,nc+2);
%  display_D_elem(gr,'colorbar',D1,[],[],'vert',8)
  minx=min(D1.dat(:));
  maxx=max(D1.dat(:));
  wx=maxx-minx;
  minx=minx-wx*0.05;
  maxx=maxx+wx*0.05;
  gray=0.1;
  for k=1:nc
    subplotgrid(gr,3:5,k+1);
    axis([ minx maxx 0.5 size(D1.dat,1)+0.5]);
    [ss,si]=sort(posterior(k,:));
    for i=1:size(D1.dat,2)
      s=si(i);
      p=posterior(k,s);
      c=[ (1-gray)*(1-p)  (1-gray) (1-gray)*(1-p)];
      if pred(s)~=tr(s) 
        if pred(s)==k
          c=c([2 1 3]);
        elseif tr(s)==k
          c=c([1 3 2]);          
        end
      end
      lh(s,k)=line(D1.dat(:,s),1:size(D1.dat,1),'Color',c);
    end
    x=[cls.mu{k}-cls.sig{k} cls.mu{k}+cls.sig{k}];
    x=(x-repmat(D1.preproc_center,1,2))./repmat(D1.preproc_scale,1,2);
    for i=1:size(D1.dat,1)
      line(x(i,:),[i i],'Color','black','LineWidth',1);
      hold on;
      plot(mean(x(i,:)),i,'k.','MarkerSize',6);
    end
    axis on
    axis ij
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    box on
  end  
end

