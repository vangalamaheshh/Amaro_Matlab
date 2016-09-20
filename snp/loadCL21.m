nms={'dat','chr','pos','sdesc','chrn','sis','supacc','supdesc','supdat','used_normals','cbs_rl','raw','medians','marker','gorigidx'};
CL21=[];
for i=1:length(nms)
  if exist(['CL21.' nms{i} '.mat'],'file')
    load(['CL21.' nms{i} '.mat']);
    eval(['CL21.' nms{i} '=' nms{i} ';']);
    eval(['clear ' nms{i} ]);
    disp(['read ' nms{i}]);
  end
end
return


load CL21.dat.mat 
 CL21.dat=dat;
%load CL21.affy_calls.mat
%load CL21.marker.mat 
load CL21.chr.mat  
load CL21.pos.mat 
load CL21.sdesc.mat 
load CL21.chrn.mat 
load CL21.history.mat 
load CL21.gorigidx.mat
if exist('CL21.sis.mat', 'file')
load CL21.sis.mat
load CL21.supacc.mat 
load CL21.supdesc.mat
load CL21.supdat.mat 
load CL21.used_normals.mat 
load CL21.temp.mat 
load CL21.cbs.mat 
load CL21.cbs_rl.mat 
load CL21.cbs_fixed.mat 
end

load CL21.raw.mat 
%load CL21.sm1.mat 
%load CL21.sm2.mat 
%load CL21.sm2j.mat 
%load CL21.sm3.mat 
%load CL21.smooth.mat 
load CL21.medians.mat

 %CL21.affy_calls=affy_calls;
% CL21.marker=marker;
 CL21.chr=chr;
 CL21.pos=pos;
 CL21.sdesc=sdesc;
 CL21.chrn=chrn;
 CL21.history=history;
 CL21.gorigidx=gorigidx;
 if exist ('sis', 'var')
     CL21.sis=sis;

 CL21.supacc=supacc;
 CL21.supdesc=supdesc;
 CL21.supdat=supdat;
 CL21.used_normals=used_normals;
 CL21.temp=temp;
 CL21.cbs=cbs;
 CL21.cbs_rl=cbs_rl;
 CL21.cbs_fixed=cbs_fixed;
 end
 CL21.raw=raw;
% CL21.sm1=sm1;
% CL21.sm2=sm2;
% CL21.sm2j=sm2j;
% CL21.sm3=sm3;
% CL21.smooth=smooth;
 CL21.medians=medians;

%clear dat affy_calls marker pos chr sdesc chrn history gorigidx sis supacc supdesc supdat used_normals temp cbs cbs_rl cbs_fixed raw sm1 sm2 sm2j sm3 smooth medians
clear dat affy_calls marker pos chr sdesc chrn history gorigidx sis supacc supdesc supdat used_normals temp cbs cbs_rl cbs_fixed raw medians
