function p=find_peaks(c,th)
%FIND_PEAKS locate peaks of data vector C.
%
%   P = FIND_PEAKS(C,TH) finds the peaks of vector C.  TH is the (required)
%   threshold for peak detection.  A peak is detected at locations where
%   the first difference is less than TH.  (i.e. raising the threshold
%   makes peaks harder to detect)  P is a vector giving the indices of the
%   peaks in C. 
%
%           Revisions
%               11 Oct 07 -- Help and documentation. -- Jen Dobson
%               (jdobson@broad.mit.edu)
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

if size(c,1)>1
  c=c';
end

dc=diff(c);
dc(abs(dc)<=th)=0;
dc(dc>th)=1;
dc(dc<-th)=-1;

rl=runlength(dc);  %find regions of 1 (+slope), 0(0 slope), -1(- slope)

if size(rl,1)==1  %c is monotonic => no peaks
  p=[];
  return;
end

if rl(1,3)==0  %don't care if first segment of c is flat
  rl=rl(2:end,:);
end

if rl(end,3)==0  %don't care if last segment of c is flat
  rl=rl(1:(end-1),:);
end

rl=rl(find(rl(:,3)~=0),:);   % remove flat parts
pk=find(rl(:,3)==1);     % segments with positive slope
pk=setdiff(pk,size(rl,1));    % reject if at end
pk2=find(rl(pk+1,3)==-1);   %
pk=pk(pk2); % the up-slopes that are followed by down-slopes (i.e. PEAKS!)


p=zeros(length(pk),1);
for i=1:length(pk)  
  range=rl(pk(i),2):rl(pk(i)+1,1);  % indices of ith peak
  [dum,pi]=min(abs(diff(c(range))));    % where does the peak peak?
  p(i)=range(pi);  %index of ith peak
end


