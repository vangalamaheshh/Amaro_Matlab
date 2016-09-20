load C2.chr.mat
load C2.chrn.mat
if exist('C2.dat1.mat','file')
  load C2.dat1.mat
  load C2.dat2.mat
else
  load C2.dat.mat
end
load C2.gorigidx.mat
load C2.history.mat
load C2.marker.mat
load C2.pos.mat
load C2.sdesc.mat
load C2.sis.mat
load C2.supacc.mat
load C2.supdat.mat
load C2.supdesc.mat
load C2.temp.mat
load C2.used_normals.mat


 C2.chr = chr;
 C2.chrn = chrn;
 if exist('C2.dat1.mat','file')
   C2.dat(1:size(dat1,1),:) = dat1;
   C2.dat(size(dat1,1)+(1:size(dat2,1)),:) =dat2;
 else
   C2.dat = dat;
 end
 C2.gorigidx = gorigidx;
 C2.history = history;
 C2.marker = marker;
 C2.pos = pos;
 C2.sdesc = sdesc;
 C2.sis = sis;
 C2.supacc = supacc;
 C2.supdat = supdat;
 C2.supdesc = supdesc;
 C2.temp = temp;
 C2.used_normals = used_normals;
 
 
 clear chr chrn dat gorigidx history marker pos sdesc sis supacc supdat supdesc temp used_normals
 
%  CN = C2;
% rewrite_data=0;
% run_on_fly=1;
% lsf_queue= 'normal';
% wait_for_jobs=1;

