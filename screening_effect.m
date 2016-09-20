function res=screening_effect(prefix,n,fgen,fscore,lognstd,T1,T2,use_mean)

prefix=[ prefix '_' num2str(lognstd)];

for ti=1:1 % T1 % FIX THIS!!!!
  factor=lognrnd(0,log(lognstd),size(n,1),1);
%  hist(log(factor),50);
  if use_mean
    factor=factor./mean(factor);
  end
  
  Fgen=factor*fgen;
  Fscore=repmat(fscore,size(n,1),1);

  K=11;
  [p,vp,x,tab]=sim_step(n*K,Fgen,Fscore,[],T2,0);
  % save screening.mat factor x
  
  g1a=squeeze(any(x,2));
  tmp=repmat(factor,1,size(g1a,2));
  tmp(g1a==0)=NaN;
  if use_mean
    f1=nanmean(tmp);
  else
    f1=nanmedian(tmp);
  end
  % typical f
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
  K1=24;
  [p1,vp1,x1,tab1]=sim_step(n1*K1,Fgen(g1,:),Fscore(g1,:),[],T2);

  g2a=squeeze(any(x1,2));
  tmp=repmat(factor(g1),1,size(g2a,2));
  tmp(g2a==0)=NaN;
  if use_mean
    f2=nanmean(tmp);
  else
    f2=nanmedian(tmp);  
  end
% typical f
  [f2t,f2i]=min(abs(f2-median(f2)));
  g2=g1(find(g2a(:,f2i)));

  tmp=repmat(n(g1,end),1,size(g2a,2));
  tmp(g2a==0)=NaN;
  if use_mean
    sz2=nanmedian(tmp);
  else
    sz2=nanmean(tmp);
  end
end % ti

save([ prefix '_screening.mat'],'factor','x','x1', 'g1', 'g2', 'f1', 'f2','sz1','sz2');


if (0)
  figure(1); clf;
  comp_dist({factor,factor(find(g1a(:,fi))),factor(g1(find(g2a(:,f2i))))},50,1)
  title(['95% CI factor(' num2str(lognstd) ')=[' ...
         num2str([exp(prctile(log(f2),2.5)),exp(prctile(log(f2),97.5))]) '], med=' num2str(exp(median(log(f2)))) ]);
  print_D([ prefix '_screening' ],{{'epsc'},{'pdf'}});
end
res={[prctile(factor,[10 25 50 75 90]);...
      prctile(factor(g1),[10 25 50 75 90]);...
      prctile(factor(g2),[10 25 50 75 90])], median(factor), median(f1), median(f2), factor, f1 ,f2, sz1, sz2};

