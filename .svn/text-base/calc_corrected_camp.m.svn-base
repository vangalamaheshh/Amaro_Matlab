function res=calc_corrected_camp(bf,V,exact_flag,tib_method,screen_flag,MYX,llrt_flag)

if ~exist('exact_flag','var')
  step=0.001; % 0.01;
  maxval=25; % was 25
  epsval=1e-20;
  exact_flag=0;
elseif isstruct(exact_flag)
  step=exact_flag.step;
  maxval=exact_flag.maxval;
  epsval=exact_flag.epsval;
  exact_flag=0;
end

if ~exist('tib_method','var')
  tib_method=0;
end
if ~exist('screen_flag','var')
  screen_flag=0;
end

if ~exist('llrt_flag','var')
  llrt_flag=0;
end

% bf={f_breast,f_colon};
[ bf{2}(1) 7.7272e-06]

for k=1:2
  p=V{k}.dat(:,1:7)./repmat(sum(V{k}.dat(:,1:7),2),1,7);
  p=[p(:,1) p(:,2)+p(:,3) p(:,4:end) ones(size(p,1),1)];
  if ~screen_flag
    vn=round(p.*repmat(sum(V{k}.dat(:,17:18),2),1,size(p,2)));
  else
    vn=round([p p].*[repmat(V{k}.dat(:,17),1,size(p,2)) repmat(V{k}.dat(:,18),1,size(p,2))]);
  end
  x=V{k}.dat(:,9:15);
  if exist('MYX','var') && ~isempty(MYX)
    x=MYX{k}(:,1:7)+MYX{k}(:,8:14);
  end
  if tib_method
    X=sum(x,2);
    pv=[];
    for j=1:size(X,1)
      if tib_method==2
        pv(j)=1-poisscdf(X(j)-1,vn(j,:)*bf{k}');
      else
        tmp=1-poisscdf(X(j)-1,vn(j,:)*bf{k}');
        pv(j)=1-sum(binoconvpdf(0:(X(j)-1),vn(j,:),bf{k}));
      end
    end
    sc=[];
    vcamp=[];
    tots=[];
  else    
    bfk=repmat(bf{k},size(vn,1),1);
    if ~screen_flag
      if ~llrt_flag
        b=ln_binopdf(x,vn,bfk);
      else
        b0=ln_binopdf(x,vn,bfk);
        b1=ln_binopdf(x,vn,x./vn);
        b1(isnan(b1))=0;
        b=b0-b1;
      end
    else
      b=ln_binopdf(MYX{k},vn,[bfk bfk]);
    end
    sc=-sum(b,2);
    vcamp=-log10(exp(-sc)*13023./(1:length(sc))');
    pv=[];
    NK=12;
    tots=nan(size(vn,1),NK+1);
    for j=1:size(vn,1); 
      disp([ j sc(j)]);
      if exact_flag
        tot=0;
        ncomb=zeros(1,NK+1);
        for kk=0:NK
          if ~screen_flag
            y=sum_combin_global(kk,7);
            if ~llrt_flag
              b=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),repmat(bf{k},size(y,1),1));
              scb=-sum(b,2);
              p=exp(-scb);
            else
              b0=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),repmat(bf{k},size(y,1),1));
              b1=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),y./repmat(vn(j,:),size(y,1),1));
              b1(isnan(b1))=0;
              b=b0-b1;
              scb=-sum(b,2);
              p=exp(sum(b0,2));
            end
            tots(j,kk+1)=sum(p(find(scb<sc(j))));
            tot=tot+tots(j,kk+1);
            ncomb(kk+1)=length(find(scb<sc(j)));
            disp([ V{k}.gacc{j} ': adding ' num2str(ncomb(kk+1)) ' comb. w/ prob ' ...
                   num2str(tots(j,kk+1)) ...
                   '=' num2str(sum(tots(j,1:kk+1))) ' pv=' num2str(1-sum(tots(j,1:kk+1)))]);
            if ncomb(kk+1)==0
              break;
            end
          else
            if (0)
              y=sum_combin_global(kk,14);
              %              b=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),repmat([bf{k} bf{k}],size(y,1),1));
              b=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),repmat([bf{k} bf{k}],size(y,1),1));
              scb=-sum(b,2);
              p=exp(-scb);
              scb(find(sum(y(:,1:7),2)==0))=0;
              scb(find(sum(y(:,8:14),2)==0))=0;
              tots(j,kk+1)=sum(p(find(scb<sc(j))));
              tot=tot+tots(j,kk+1);
              ncomb(kk+1)=length(find(scb<sc(j)));
              from_zero_scb=sum(p(find(scb==0)));
              disp([ 'adding ' num2str(ncomb(kk+1)) ' comb. w/ prob ' ...
                     num2str(tots(j,kk+1)) ...
                     ' (' num2str(from_zero_scb) ',' num2str(tots(j,kk+1)-from_zero_scb)  ')' ...
                     '=' num2str(sum(tots(j,1:kk+1)))]);
              if tots(j,kk+1)<1e-10
                disp(['Stopped after n_mut=' num2str(kk) ' : ' num2str(tots(j,kk+1))]);
                break
              end
            else
              if kk==0
                p_zero_in_D=exp(vn(j,1:7)*(log(1-bf{k}))');
                p_zero_in_V=exp(vn(j,8:14)*(log(1-bf{k}))');
                tots(j,1)=p_zero_in_D+p_zero_in_V-p_zero_in_D*p_zero_in_V;
                disp([p_zero_in_D p_zero_in_V tots(j,1)]);
              elseif kk<2
                tots(j,kk+1)=0;
              else
                if (0)
                  y_D=sum_combin_global(1,7);
                  y_D=[ y_D zeros(size(y_D,1),7)];
                  y_V=sum_combin_global(1,7);
                  y_V=[ zeros(size(y_V,1),7) y_V];              
                  y=sum_combin_global(kk-2,14);
                  c=counter([ size(y,1) size(y_D,1) size(y_V,1)]);
                  y=y(c(:,1),:)+y_D(c(:,2),:)+y_V(c(:,3),:);
                  y=unique(y,'rows');
                else
                  y=sum_combin_global(kk,14);
                  y(find(sum(y(:,1:7),2)==0),:)=[];
                  y(find(sum(y(:,8:14),2)==0),:)=[];
                end
                %              b=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),repmat([bf{k} bf{k}],size(y,1),1));
                b=ln_binopdf(y,repmat(vn(j,:),size(y,1),1),repmat([bf{k} bf{k}],size(y,1),1));
                scb=-sum(b,2);
                p=exp(-scb);
                tots(j,kk+1)=sum(p(find(scb<sc(j))));
                ncomb(kk+1)=length(find(scb<sc(j)));
                disp([ 'adding ' num2str(ncomb(kk+1)) ' comb. w/ prob ' ...
                       num2str(tots(j,kk+1)) '=' num2str(sum(tots(j,1:kk+1)))]);
              end
              tot=tot+tots(j,kk+1);
              if kk>=2 && tots(j,kk+1)<1e-10
                disp(['Stopped after n_mut=' num2str(kk) ' : ' num2str(tots(j,kk+1))]);
                break
              end
            end
          end
        end
        if ncomb(end)>0
          disp(['WARNING: ' num2str(k) ':' num2str(j) ' - ' num2str(tots(j,kk)) ' ' num2str(tots(j,kk+1))]);
        end
        pv(j)=1-tot;
      else
        d=exact_dist_pdf(vn(j,:),bf{k},0:step:maxval,epsval,50); %min(maxval,sc(j)+step)
        if sc(j)>maxval
          pv(j)=0;
        else
%          pv(j)=sum(d((floor(sc(j)/step)+1):end));
          bins=floor(-b(j,:)/step)+1;
          curbin=sum(bins)-length(bins)+1;
%          tmp=sum(d((floor(sc(j)/step)+1):end));
%          pv(j)=1-sum(d(1:floor(sc(j)/step))); 
          tmp=sum(d(curbin:end));
          pv(j)=1-sum(d(1:(max(0,curbin-1))));
          disp([ 'diff: ' num2str(pv(j)-tmp)]);
        end
      end
    end
  end
  [spv,spi]=sort(pv);
  qv=spv*13023./(1:length(pv));
  for xx=(length(qv)-1):-1:1
    qv(xx)=min(qv(xx),qv(xx+1));
  end
  camp=-log10(qv);
  res{k}={bf,camp,pv,spi,qv,sc,vcamp,tots};
end

