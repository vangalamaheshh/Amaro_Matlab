function [idx,res]=outliers(x,params)

if ischar(params)
  tmp.method=params;
  params=tmp;
end
res=[];

switch params.method
 case 'medmad',
  med=median(x);
  s=mad(x,1)*(1/norminv(0.75));
  res.rng=[med+params.mads*s med-params.mads*s];
  idx=find(x>res.rng(2) | x<res.rng(1));
 case 'iqr' 
  % param.m = number of iqr's beyond which to find outliers
  if ~isfield(params,'n'), params.n=1; end;
  pc=prctile(x,[25 75]);
  iqr=diff(pc);
  res.rng=[pc(1)-(params.n*iqr) pc(2)+(params.n*iqr)];
  idx=find(x>res.rng(2) | x<res.rng(1));
 case 'copa'
  % param.p = percentile cutoff
  x=x-repmat(median(x,2),1,size(x,2));
  x=x./repmat(mad(x,1,2),1,size(x,2));
  score=prctile(x,params.p,2);

  [s,idx]=sort(score,1,'descend');
  res.score=score;
end
