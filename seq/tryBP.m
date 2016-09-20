function X = tryBP(individual,X,idx)
% tryBP(individual,X[,idx])
%
% real-time interface to BreakPointer (Yotam Drier, 2010)
%
% individual = e.g. GBM-0188
% X = a dRanger results struct
% idx = (optional) which records in X to try in BP
%
% Mike Lawrence 2010-02-10

while true
  rand('twister',cputime);
  tempdir = ['/tmp/' num2str(round(1e9*rand))];
  if ~exist(tempdir,'dir'), break; end
end

fprintf('Running in %s\n',tempdir);
tdir = [tempdir '/tumor'];
ndir = [tempdir '/normal'];
fhname = ['/xchip/cga1/firehose_output/Individual/' individual '/wgs'];
tsamp = [individual '_Tumor'];
nsamp = [individual '_Normal'];
tbam = [fhname '/bam/tumor.bam'];
nbam = [fhname '/bam/normal.bam'];
blacklist = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/reference/lane_blacklist.txt';
refdir = '/xchip/tcga/gbm/analysis/lawrence/genome/hg18';
libdir = '/xchip/tcga/gdac_prod/applications/genepattern/taskLib/BreakPointer.41.1321/';
gscript = ['/home/radon00/lawrence/CancerGenomeAnalysis/trunk/analysis_pipeline/genepattern/modules/BreakPointer/' ...
           'breakpointer_gather.pl'];
keeptempfiles = 1;

p = pwd;
mkdir(tempdir);
drr = [tempdir '/drr.txt'];
if exist('idx','var'), X = reorder_struct(X,idx); end
save_struct(X,drr);
mkdir(tdir);cd(tdir);
tsc = [tdir '/scatter.1'];
mkdir(tsc);cd(tsc);
fh_BreakPointerScatter(tsamp,drr,tbam,blacklist,refdir,libdir,1,slength(X));
cd('..');
system(['perl ' gscript ' ' tsamp ' ' libdir ' . ' num2str(keeptempfiles)]);
mkdir(ndir);cd(ndir);
nsc = [ndir '/scatter.1'];
mkdir(nsc);cd(nsc);
fh_BreakPointerScatter(nsamp,drr,nbam,blacklist,refdir,libdir,1,slength(X));
cd('..');
system(['perl ' gscript ' ' nsamp ' ' libdir ' . ' num2str(keeptempfiles)]);
cd(p);

BPt = [tdir '/' tsamp '.breakpoints.txt'];
BPn = [ndir '/' nsamp '.breakpoints.txt'];

X = make_numeric(X,{'chr1','chr2','pos1','pos2','str1','str2'});

flds = {'line','num','chr1','BPpos1','chr2','BPpos2','inversion','SWreads','SWscore',...
  'firstseq','lenhomology','lenforeign','foreignseq','BWAreads'};
flds2 = setdiff(flds,{'line','num','foreignseq'});
T = make_numeric(rename_fields(load_struct_noheader(BPt),colx(1:14),flds),flds2);
N = make_numeric(rename_fields(load_struct_noheader(BPn),colx(1:14),flds),flds2);
if isnumeric(X.num), num = num2cellstr(X.num); else num = X.num; end
tidx = listmap(T.num,num);
nidx = listmap(N.num,num);
if any(T.chr1~=X.chr1(tidx)) || any(T.chr2~=X.chr2(tidx)) || any (T.inversion~=(X.str1(tidx)==X.str2(tidx)))
  error('%s does not match X\n',BPt);
end
if any(N.chr1~=X.chr1(nidx)) || any(N.chr2~=X.chr2(nidx)) || any (N.inversion~=(X.str1(nidx)==X.str2(nidx)))
  error('%s does not match X\n',BPn);
end
T.BPhit = ones(slength(T),1);
N.BPhit = ones(slength(N),1);
T.diffpos1 = nan(slength(T),1);  T.diffpos2 = nan(slength(T),1);
N.diffpos1 = nan(slength(N),1);  N.diffpos2 = nan(slength(N),1);
flds3 = {'BPhit','BPpos1','diffpos1','BPpos2','diffpos2','SWreads','SWscore',...
  'firstseq','lenhomology','lenforeign','foreignseq','BWAreads'};
T = keep_fields(T,flds3);  N = keep_fields(N,flds3);
X = struct_assign(X,tidx,rename_fields(T,flds3,regexprep(flds3,'(.*)','T_$1')));
X = struct_assign(X,nidx,rename_fields(N,flds3,regexprep(flds3,'(.*)','N_$1')));
X.T_BPhit = (X.T_BPhit==1);
X.N_BPhit = (X.N_BPhit==1);
X.T_diffpos1 = X.T_BPpos1-X.pos1;  X.T_diffpos2 = X.T_BPpos2-X.pos2;
X.N_diffpos1 = X.N_BPpos1-X.pos1;  X.N_diffpos2 = X.N_BPpos2-X.pos2;

