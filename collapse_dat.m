function [cdat,ni]=collapse_dat(dat,vec,vecval,method)
% [cdat,ni]=collapse_dat(dat,vec,vecval,method)
%    collapses the rows of dat according to the values of vec
%    the output cdat has rows as in vecval. Each row is 
%    the aggregation of all rows in dat whose vec value is that of 
%    the corresponding row in vecval.
%    method - indicates the method to perform the aggregation: 
%             sum, max, any, mean
%
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

if ischar(method)
  tmp.method=method;
  method=tmp;
end

ni=zeros(length(vecval),1);
cdat=NaN*ones(length(vecval),size(dat,2));
idx=1:length(vecval);
if (1) % lines with a unique val can be copied over in many methods
  switch method.method
   case {'mean','nanmean','sum','max','min','median','nanmedian','mode'}
    [u,ui,uj]=unique(vec);
    hc=histc(uj,1:length(ui));
    one_idx=find(ismember(uj,find(hc==1)));
    [Mt,m1,m2]=match_num_sets(vec(one_idx),vecval);
    cdat(m2,:)=dat(one_idx(m1),:);
%    disp(length(m2));
    idx=setdiff(idx,m2);
  end
end

for i=idx % 1:length(vecval)
  if mod(i,1000)==0
    disp(i);
  end
  pos=find(vec==vecval(i));
  ni(i)=length(pos);
  if ni(i)>0
    switch method.method
     case 'sum'
      cdat(i,:)=sum(dat(pos,:),1);
     case 'max'
      cdat(i,:)=max(dat(pos,:),[],1);
     case 'any'
      cdat(i,:)=any(dat(pos,:),1);
     case 'mean'
      cdat(i,:)=mean(dat(pos,:),1);
     case 'nanmean'
      cdat(i,:)=nanmean(dat(pos,:),1);
     case 'median'
      cdat(i,:)=median(dat(pos,:),1);
     case 'nanmedian'
      cdat(i,:)=nanmedian(dat(pos,:),1);
     case 'mode'
      cdat(i,:)=mode(dat(pos,:),1);      
     case 'maxscore'
      [mx,mi]=max(method.score(pos));
      cdat(i,:)=dat(pos(mi),:);
     case 'minscore'
      [mn,mi]=min(method.score(pos));
      cdat(i,:)=dat(pos(mi),:);
    end
  end
end

    

