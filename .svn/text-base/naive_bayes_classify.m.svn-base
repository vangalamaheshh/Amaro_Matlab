function [class,score,sig,logscore]=naive_bayes_classify(pts,trpts,cls,w,tr)
% [class,score,sig,logscore]=naive_bayes_classify(pts,mempts,cls,tr)
%    performs a naive-Bayes classifier.
%    pts - data points to be classified (rows)
%    trpts - training data points (rows)
%    cls - class index of training set (matrix n_cls x n_pts) 
%    w - weight of each class
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

% train model

nc=size(cls,1);
ntr=size(cls,2);
npts=size(pts,1);
c=cls./repmat(sum(cls,2),1,ntr);
for c=1:size(cls,1)
  models.mu=c*trpts;
  s2=c*(trpts.^2);
  models.inv_sig2=1./((s2-model(c))./(ntr-1)).^2;
end

logscore=zeros(1,size(pts,1));
expon=zeros(npts,nc);
for c=1:nc
  %%%%%%%%%%%%%%%%
  d2=dist(models.mu(c,:),pts,struct('method','seuclidean_sqrd','',1./(models.sig(c,:)).^2*ones(1,3)));
  p_x_giv_c=exp(expon(:,c)-(d2/2);
end

for i=1:size(pts,1)
  end
  expon=-(d.^2/(2*sig.^2));
    expon=expon-expon(1);
    f=exp(expon);
    c=cls(:,di);
    % t=sum(c,2); % t(2)/t(1)
    c=c./repmat(sum(c,2),1,size(c,2));
    ck=(c*f).*w;
    ck=ck/sum(ck,1);
    [score(i),class(i)]=max(ck);
    if exist('tr','var')
      logscore(i)=tr(:,i)'*log((ck+eps)/(1+eps));
    end
end



