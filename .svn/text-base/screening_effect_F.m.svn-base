function RES=screening_effect_F(prefix,n,fgen,fscore,lognstd,T1,T2,use_mean,K,K1,notab)


if ~exist('K','var')
  K=11;
end
if ~exist('K1','var')
  K1=24;
end

if ~exist('notab','var')
  notab=0;
end

RES=cell(T1,14);
for ti=1:T1 
  disp(['ITERATION=' num2str(ti)]);
  prefix=[ prefix '_' num2str(ti) '_' num2str(lognstd)];
  factor=lognrnd(0,log(lognstd),size(n,1),1);
%  hist(log(factor),50);
  if use_mean
    factor=factor./mean(factor);
  end
  
  Fgen=factor*fgen;
  Fscore=repmat(fscore,size(n,1),1);

  [p,vp,x,tab]=sim_step(n*K,Fgen,Fscore,[],T2,0,notab);
  % save screening.mat factor x
  
  g1a=squeeze(any(x,2));
  nmut1=squeeze(sum(sum(x,1),2));
  ngene1=sum(g1a,1);
  tmp=repmat(factor,1,size(g1a,2));
  tmp(g1a==0)=NaN;
  if use_mean
    f1=nanmean(tmp);
  else
    f1=nanmedian(tmp);
  end
  % find typical factor, ft
  [ft,fi]=min(abs(f1-median(f1)));
  
  % size effect
  tmp=repmat(n(:,end),1,size(g1a,2));
  tmp(g1a==0)=NaN;
  if use_mean
    sz1=nanmean(tmp);
  else
    sz1=nanmedian(tmp);
  end
%comp_dist({factor,factor(find(g1a(:,fi)))},50,1)
%title(['90% CI factor(' num2str(lognstd) ')=[' ...
%       num2str([exp(prctile(log(f1),10)),exp(prctile(log(f1),90))]) '], med=' num2str(exp(median(log(f1)))) ]);
  g1=find(g1a(:,fi));

  n1=n(g1,:);
  [p1,vp1,x1,tab1]=sim_step(n1*K1,Fgen(g1,:),Fscore(g1,:),[],T2,0,notab);

  g2a=squeeze(any(x1,2));
  nmut2=squeeze(sum(sum(x1,1),2));
  ngene2=sum(g2a,1);
  tmp=repmat(factor(g1),1,size(g2a,2));
  tmp(g2a==0)=NaN;
  if use_mean
    f2=nanmean(tmp);
  else
    f2=nanmedian(tmp);  
  end
% find typical f after second round, f2t
  [f2t,f2i]=min(abs(f2-median(f2)));
  g2=g1(find(g2a(:,f2i)));

  tmp=repmat(n(g1,end),1,size(g2a,2));
  tmp(g2a==0)=NaN;
  if use_mean
    sz2=nanmean(tmp);
  else
    sz2=nanmedian(tmp);
  end
%  RES(ti,:)={factor,sparse(squeeze(sum(x,2))),sparse(squeeze(sum(x1,2))),g1,g2,f1,f2,sz1,sz2};
  xx=x(g1,:)+x1;
  g2=find(any(x1,2));
  VT=binopdf(xx(g2,:),n(g1(g2),:)*35,Fscore(g1(g2),:));
  vp=prod(VT,2);
  [svp,svpi]=sort(vp);
  vcamp=nan(length(vp),1);
  vcamp(svpi)=-log10(svp*13023./(1:length(svp))');
  cang=g1(g2(find(vcamp>1)));
  ncan=length(cang);
  if use_mean
    szg=nanmean(n(cang,end));
  else
    szg=nanmedian(n(cang,end));
  end
  est_rates=sum(x)./sum(n)/K;
  RES(ti,:)={nmut1,ngene1,nmut2,ngene2,g1,g2,f1,f2,sz1,sz2,ncan,cang,szg,est_rates};
%   save([ prefix '_screening.mat'],'factor','x','x1', 'g1', 'g2', 'f1', 'f2','sz1','sz2');
end % ti



if (0)
  figure(1); clf;
  comp_dist({factor,factor(find(g1a(:,fi))),factor(g1(find(g2a(:,f2i))))},50,1)
  title(['95% CI factor(' num2str(lognstd) ')=[' ...
         num2str([exp(prctile(log(f2),2.5)),exp(prctile(log(f2),97.5))]) '], med=' num2str(exp(median(log(f2)))) ]);
  print_D([ prefix '_screening' ],{{'epsc'},{'pdf'}});
end

if (0)
  res={[prctile(factor,[10 25 50 75 90]);...
        prctile(factor(g1),[10 25 50 75 90]);...
        prctile(factor(g2),[10 25 50 75 90])], median(factor), median(f1), median(f2), factor, f1 ,f2, sz1, sz2, g1, ...
       g2};
end  
