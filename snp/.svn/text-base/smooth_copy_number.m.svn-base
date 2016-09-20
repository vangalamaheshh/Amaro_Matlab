function [C,res]=smooth_copy_number(C,sz,padval,sm_method)
%SMOOTH_COPY_NUMBER 
%
%    [C,RES] = SMOOTH_COPY_NUMBER(C,SZ,PADVAL,SM_METHOD)
%
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

if ischar(sm_method)
  sm_method1.method=sm_method;
  sm_method=sm_method1;
end

if ~isfield(C,'chrn')
  C.chrn=chromosome2num(C.chr);
end
res=[];
u_chr=unique(C.chrn);
if strcmp(class(C),'datastruct') && isdiskfield(C,'dat')
    filename = regexprep(get_datafile(C,'dat'),'dat','smooth');
C = add_diskfield(C,filename,'smooth',getsize(C),'single');
end


for i=1:length(u_chr)
  f=find(C.chrn==u_chr(i));
  switch sm_method.method
   case {'mean','median','min','prodrange'}
    C.smooth(f,:)=running_window(C.dat(f,:),sz,padval,sm_method);
   case 'segmented_mean'
    for k=1:size(C.dat,2)
      sm_method.segments=segment_copy_number(C.dat(f,k),sm_method.seg_sz, ...
                              sm_method.diff_method, ...
                              sm_method.pcutoff);
      seg(f,k)=sm_method.segments';
%      keyboard
      C.smooth(f,k)=running_window(C.dat(f,k),sz,padval,sm_method);
      k
    end
    res=seg;
   otherwise
    error('no such method');
  end
  verbose(['finished chromosome: ' num2str(u_chr(i))],10);
end


