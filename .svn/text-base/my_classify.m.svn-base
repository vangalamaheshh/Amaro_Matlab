function [confmat,nerr]=my_classify(D,method,tps)

D.supdat=D.supdat(tps,:);
D.supdat(isnan(D.supdat))=0;
if length(find(sum(D.supdat,1)==0))>0 
  D.supdat=[D.supdat; (sum(D.supdat,1)==0)];
end
ntp=size(D.supdat,1);

n=size(D.dat,2);

confmat=zeros(n,ntp,ntp);
nerr=zeros(n,1);
for i=1:n
  conf=zeros(ntp,ntp);
  looD=D;
  trn_set=setdiff(1:n,i);
  tst_set=i;
  looD.dat=looD.dat(:,trn_set);
  looD.sdesc=looD.sdesc(trn_set,:);
  looD.supdat=looD.supdat(:,trn_set);
  looD.affy_call=looD.affy_call(:,trn_set);
  
  ntrain=length(trn_set);
  prior=sum(looD.supdat,2)+ones(size(looD.supdat,1),1);
  
  if isfield(method,'weigh_features')
    fet_w=weigh_features(looD,method.weigh_features,ones(size(prior)));
  end
  switch method.type
      case 'knn'
          [ind,tmp]=find(looD.supdat);
          class=knn_classify(D.dat(:,tst_set)',D.dat(:,trn_set)', ...
                             looD.supdat,fet_w,ones(size(prior)),method.k); % 1./prior
      case 'naive'
      case 'classify'
          [ind,tmp]=find(looD.supdat);
          keyboard
          class=classify(D.dat(:,tst_set)',D.dat(:,trn_set)',ind,'linear',prior'/ntrain);
  end
  true=find(D.supdat(:,i));
  if length(true)>1
      disp(['Ambiguous type at point ' num2str(i)]);
      true=true(1);
  end
  conf(class,true)=conf(class,true)+1;
  cnerr=sum(sum(conf.*(ones(ntp,ntp)-eye(ntp))));
  nerr(i)=cnerr;
  confmat(i,:,:)=conf;
end
