function [X,TT,NN] = load_dR_and_BP_data(individual,dRmatfile,BPt,BPn,diffthresh)

% load dRanger results
load(dRmatfile,'X')
if ~exist('X','var'), error('%s does not contain X',dRmatfile); end
if nargout>=2
  load(dRmatfile,'TT')
  if ~exist('TT','var'), error('%s does not contain TT',dRmatfile); end
end
if nargout>=3
  load(dRmatfile,'NN')
  if ~exist('NN','var'), error('%s does not contain NN',dRmatfile); end
end


if ~exist('diffthresh','var'), diffthresh = 100; end

% check for case of NO DATA
if isfield(X,'NO_DATA'), return; end

% copy dRanger's original positions to new fields

X.dRpos1 = X.pos1;
X.dRpos2 = X.pos2;

% load BreakPointer results

orig_flds = {'bpnum','drnum','chr1','pos1','chr2','pos2','inversion','bpreads','qual','seq2first',...
  'exact_mh','fs_len','fs','bwareads','soft_mh'};
flds = {'line','num','chr1','BPpos1','chr2','BPpos2','inversion','SWreads','SWscore',...
  'firstseq','lenhomology','lenforeign','foreignseq','BWAreads','lenhomology_soft'};
flds2 = setdiff(flds,{'line','num','foreignseq'});
T = load_struct(BPt);
N = load_struct(BPn);
T = rename_fields(T,orig_flds,flds);
N = rename_fields(N,orig_flds,flds);
T = make_numeric(T,flds2);
N = make_numeric(N,flds2);

if isnumeric(X.num), num = num2cellstr(X.num); else num = X.num; end
tidx = listmap(T.num,num);
nidx = listmap(N.num,num);
if any(T.chr1~=X.chr1(tidx)) || any(T.chr2~=X.chr2(tidx)) || any(T.inversion~=-1 & T.inversion~=(X.str1(tidx)==X.str2(tidx)))
  error('%s does not match %s\n',BPt,dRmatfile);
end
if any(N.chr1~=X.chr1(nidx)) || any(N.chr2~=X.chr2(nidx)) || any(N.inversion~=-1 & N.inversion~=(X.str1(nidx)==X.str2(nidx)))
  error('%s does not match %s\n',BPn,dRmatfile);
end
T.BPhit = ones(slength(T),1);
N.BPhit = ones(slength(N),1);
T.diffpos1 = nan(slength(T),1);  T.diffpos2 = nan(slength(T),1);
N.diffpos1 = nan(slength(N),1);  N.diffpos2 = nan(slength(N),1);
flds3 = {'BPhit','BPpos1','diffpos1','BPpos2','diffpos2','SWreads','SWscore',...
  'firstseq','lenhomology','lenhomology_soft','lenforeign','foreignseq','BWAreads'};
T = keep_fields(T,flds3);  N = keep_fields(N,flds3);
X = struct_assign(X,tidx,rename_fields(T,flds3,regexprep(flds3,'(.*)','T_$1')));
X = struct_assign(X,nidx,rename_fields(N,flds3,regexprep(flds3,'(.*)','N_$1')));
X.T_BPhit = (X.T_BPhit==1);
X.N_BPhit = (X.N_BPhit==1);
X.T_diffpos1 = X.T_BPpos1-X.pos1;  X.T_diffpos2 = X.T_BPpos2-X.pos2;
X.N_diffpos1 = X.N_BPpos1-X.pos1;  X.N_diffpos2 = X.N_BPpos2-X.pos2;

% summarize BreakPointer results
% NaN = not tried;  -1 = succeeded in normal; 0 = failed in both; 1 = succeeded in tumor, failed in normal


minswreads1 = 2;
minswreads2 = 3;
minswscore1 = 0.85;
minswscore2 = 0.8;
X.BPresult = nan(slength(X),1);
X.BPsomaticratio = nan(slength(X),1);
X.approxflag = nan(slength(X),1);
for i=1:slength(X)
  if X.BPtry(i)
    t = (X.T_BPhit(i) && ((X.T_SWreads(i)>=minswreads1 && X.T_SWscore(i)>=minswscore1)...
                          || (X.T_SWreads(i)>=minswreads2 && X.T_SWscore(i)>=minswscore2)) ...
                          && abs(X.T_diffpos1(i))<=diffthresh && abs(X.T_diffpos2(i))<=diffthresh);
    n = (X.N_BPhit(i) && ((X.N_SWreads(i)>=minswreads1 && X.N_SWscore(i)>=minswscore1)...
                          || (X.N_SWreads(i)>=minswreads2 && X.N_SWscore(i)>=minswscore2)) ...
                          && abs(X.N_diffpos1(i))<=diffthresh && abs(X.N_diffpos2(i))<=diffthresh);
    X.BPsomaticratio(i)=X.T_BWAreads(i)/X.N_BWAreads(i);    
    X.approxflag(i)=~t; 
    if n && (X.BPsomaticratio(i) < 4), X.BPresult(i) = -1;
    elseif t, X.BPresult(i) = 1;
    else X.BPresult(i) = 0; end
  end
end

% in cases where BreakPointer succeed, overwrite pos1 and pos2 with T_BPpos1 and T_BPpos2
idx = find(X.approxflag==false);
X.pos1(idx) = X.T_BPpos1(idx);
X.pos2(idx) = X.T_BPpos2(idx);
