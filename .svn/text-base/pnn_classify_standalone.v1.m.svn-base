function [class,score,sig,logscore]=pnn_classify_standalone(pts,mempts,cls,w,sig,sig_type,dist_type,tr)
% [class,score,sig,logscore]=pnn_classify(pts,mempts,cls,w,sig,dist_type,tr)
%    performs a pnn classifier.
%    pts - data points to be classified (rows)
%    mempts - training data points (rows)
%    cls - class index of training set
%    w - weight of each class
%    sig - sigma of the gaussian. 
%    sig_type - 0: use raw sigma values
%               1: use sig*(mean of positive nearest neighbor distances in the training set)
%               2: use sig*(median of positive nearest neighbor distances in the training set)
%    dist_type - distance type of Matlab's pdist function
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

% calculate sig based on a scale determined from the training set distances (sig_type > 0)
if sig_type>0
  % calculate distances among the training set
  dtrain=dist(mempts,[],dist_type); 
  % set the digonal to infinity
  dtrain=dtrain+diag(Inf*ones(size(dtrain,1),1));
  
  % find minimal distance
  min_d=min(dtrain);

  % consider only distances greater than 0 to handle cases in which the same
  % training vector is given more than once
  if sig==2
    scale=median(min_d(min_d>0)); 
  else
    scale=mean(min_d(min_d>0));    
  end

  sig=sig*scale;
end


logscore=zeros(1,size(pts,1));
% go over points to predict
for i=1:size(pts,1)
    % calculate distance to training set
    d=dist(mempts,pts(i,:),dist_type);
    [ds,di]=sort(d);
    expon=-(ds.^2/(2*sig.^2));
    expon=expon-expon(1);
    f=exp(expon);
    c=cls(:,di);
    c=c./repmat(sum(c,2),1,size(c,2));
    ck=(c*f).*w;
    ck=ck/sum(ck,1);
    [score(i),class(i)]=max(ck);
    if exist('tr','var')
      logscore(i)=tr(:,i)'*log((ck+eps)/(1+eps));
    end
end


function d=dist(X,Y,dist_type)

if isempty(Y)
  d=squareform(pdist(X,dist_type));
else
  d=dist([X; Y],[],dist_type);
  nx=size(X,1);
  d=d(1:nx,nx+(1:size(Y,1)));
end

