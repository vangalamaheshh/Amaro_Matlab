function d=running_window(dat,sz,padval,sm_method)
%RUNNING_WINDOW smooth data using a running window.
%
%   D = RUNNING_WINDOW(DAT,SZ,PADVAL,SM_METHOD) returns smoothed data D by 
%   padding data array DAT with rows of value PADVAL (usually 0).  Then 
%   data is using information from the string or structure SM_METHOD.  
%   The allowed values of SM_METHOD (or SM_METHOD.METHOD) are:
%      'segmented_mean': takes the mean over the part of the size SZ window 
%           that is within a segment.  Segments must be defined as a vector
%           with length size(dat,1) in the SM_METHOD field SM_METHOD.SEGMENTS. 
%      'mean': takes the running mean of the data, with a window of size SZ
%      'median': takes the running median of the data, window size SZ
%      'min': takes the running minimum of the data, window size SZ
%      'prodrange': the product of the running min and the running max
%      'prctile': takes the running percentile of the data.  The percentile 
%           value must be given in SM_METHOD.PRCTILE.

%% Pad and segment

hsz1=floor(sz/2);
hsz2=sz-hsz1-1;
dat=[ repmat(padval,hsz1,size(dat,2)); dat; repmat(padval,hsz2,size(dat,2))]; 
d=zeros(size(dat));

if ischar(sm_method)
  sm_method.method=sm_method;
end

switch sm_method.method
 case 'segmented_mean'
  seg=[ repmat(sm_method.segments(1),1,hsz1) sm_method.segments repmat(sm_method.segments(end),1,hsz2)];
end

%% Smooth the data

switch sm_method.method  
 case 'mean'
  for j=1:size(d,2)
    tmp=conv(dat(:,j),ones(sz,1)/sz);  %running mean
    d(:,j)=tmp((hsz2+1):(end-hsz1));
  end
 otherwise
  for i=(hsz1+1):(size(dat,1)-hsz2)  
    switch sm_method.method
     case 'median',
      d(i,:)=median(dat((i-hsz1):(i+hsz2),:),1); 
     case 'mean',
      d(i,:)=mean(dat((i-hsz1):(i+hsz2),:),1);    
     case 'min',
      d(i,:)=min(dat((i-hsz1):(i+hsz2),:),[],1);   
     case 'prctile'
      d(i,:)=prctile(dat((i-hsz1):(i+hsz2),:),sm_method.prctile);   
     case 'prodrange'
      d(i,:)=min(dat((i-hsz1):(i+hsz2),:),[],1).*max(dat((i-hsz1):(i+hsz2),:),[],1);
     case 'segmented_mean'
      window=(i-hsz1):(i+hsz2);
      d(i,:)=mean(dat(window(find(seg(window)==seg(i))),:),1);
     otherwise
      disp('no such method');
    end
  end
end

d=d(hsz1+1:(size(dat,1)-hsz2),:);
