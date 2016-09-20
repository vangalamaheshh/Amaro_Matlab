function [amp_sc,del_sc]=snp_score(C,score_method)
% [amp_sc,del_sc]=snp_score(C,score_method)
%     Calculate a score per snp (marker).
%
% ---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

if ischar(score_method)
  tmp.method=score_method;
  score_method=tmp;
end

switch score_method.method
 case 'nxa'

  if isfield(score_method,'amp_thresh')
    amp_thresh=score_method.amp_thresh;
  else  
    error('No amp_thresh');
  end
  
  if isfield(score_method,'del_thresh')
    del_thresh=score_method.del_thresh;
  else  
    error('No del_thresh');
  end
  
  x=double(C.dat);
  x(isinf(x))=NaN; % treat Inf/-Inf as NaN
  x(x<amp_thresh)=0;
%   x=sparse(x);
%  xnan=sum(isnan(x),2);
%   x(isnan(x))=0;
  
  y=double(-C.dat);
  y(isinf(y))=NaN; % treat Inf/-Inf as NaN
  y(y<-del_thresh)=0;
%  y=sparse(y);
%  ynan=sum(isnan(y),2);
%  y(isnan(y))=0;

  amp_sc=full(nanmean(x,2));
  del_sc=full(nanmean(y,2));
 
 case 'freq'

  if isfield(score_method,'amp_thresh')
    amp_thresh=score_method.amp_thresh;
  else  
    error('No amp_thresh');
  end
  
  if isfield(score_method,'del_thresh')
    del_thresh=score_method.del_thresh;
  else  
    error('No del_thresh');
  end

  x=double(C.dat);
  x(x<amp_thresh)=0;
  x(x>=amp_thresh)=1;
%  x=sparse(x);
%  xnan=sum(isnan(x),2);
%  x(isnan(x))=0;
  
  y=double(-C.dat);
  y(y<-del_thresh)=0;
  y(y>=-del_thresh)=1;
%  y=sparse(y);
%  ynan=sum(isnan(y),2);
%  y(isnan(y))=0;

  amp_sc=full(nanmean(x,2));
  del_sc=full(nanmean(y,2));
  
 otherwise
  error('No such score method');
end

