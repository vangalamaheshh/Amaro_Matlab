function [class,score,sig,logscore]=pnn_classify(pts,mempts,cls,w,sig_type,sig,dist_type,tr)
% [class,score,sig,logscore]=pnn_classify(pts,mempts,cls,w,sig,dist_type,tr)
%    performs a pnn classifier.
%    pts - data points to be classified (rows)
%    mempts - training data points (rows)
%    cls - class index of training set
%    w - weight of each class
%    sig_type - 
%          0 - use sig
%          1 - sig*(mean of positive nearest neighbor distances in the training set)
%          2 - sig*(median of positive nearest neighbor distances in the training set)
%    sig - sigma of the gaussian. 
%    dist_type - distance type
%    tr - (optional) true classification of points to calc. logscore 
%         matrix of no. of classes x no. of pts with 1 in element i,j
%         if point j belongs to class i
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%
% This software and its documentation are copyright 2005 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
% or functionality. 
% 
if sig_type>0
  dtrain=dist(mempts,[],dist_type);
  dtrain=dtrain+diag(Inf*ones(size(dtrain,1),1));
  min_d=min(dtrain);
  if sig_type==2
    scale=median(min_d(min_d>0)); 
  else
    scale=mean(min_d(min_d>0));    
  end
  sig=-sig*scale;
end

logscore=zeros(1,size(pts,1));
for i=1:size(pts,1)
    d=dist(mempts,pts(i,:),dist_type);
%    [ds,di]=sort(d);
% prev. version had:  ds=ds-ds(1);
    expon=-(d.^2/(2*sig.^2));
    expon=expon-max(expon); % expon(1);
    f=exp(expon);

    % calculate posterior probability
    c=cls./repmat(sum(cls,2),1,size(cls,2));
    ck=(c*f).*w;
    ck=ck/sum(ck,1);

    % c=cls(:,di);
    % t=sum(c,2); % t(2)/t(1)
    [score(i),class(i)]=max(ck,[],1);
    ck(ck<eps)=eps;    
    if exist('tr','var')
      logscore(i)=tr(:,i)'*log(ck); % (ck+eps)/(1+eps));
    end
end


