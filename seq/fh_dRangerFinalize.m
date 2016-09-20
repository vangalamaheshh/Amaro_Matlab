function fh_dRangerFinalize(individual,dRmatfile,circospng,BPt,BPn,diffthresh)
% fh_dRangerFinalize(individual,dRmatfile,circospng,BPt,BPn)
%
% Mike Lawrence 2010

if nargin~=6, error('Usage: fh_dRangerFinalize(individual,dRmatfile,circospng,BPt,BPn,diffthresh)'); end
if ~isnumeric(diffthresh), diffthresh=str2double(diffthresh); end

MIN_FINAL_SOMATIC_SCORE = 4;   % (need to expose this parameter to FH)

fprintf('fh_dRangerFinalize\n');
fprintf(['  individual = ' individual '\n']);
fprintf(['  dRmatfile = ' dRmatfile '\n']);
fprintf(['  circospng = ' circospng '\n']);
fprintf(['  BPt = ' BPt '\n']);
fprintf(['  BPn = ' BPn '\n']);
fprintf('  diffthresh = %d\n',diffthresh);
fprintf(['  MIN_FINAL_SOMATIC_SCORE = ' MIN_FINAL_SOMATIC_SCORE '\n']);

fprintf('Loading input data\n');
X = load_dR_and_BP_data(individual,dRmatfile,BPt,BPn,diffthresh);

% check for case of NO DATA
if isfield(X,'NO_DATA')
  H = ['<h3>dRanger results for ' individual '</h3>'...
       '<p><font color=red size=+1>NO DATA</font>'...
       '<br>Does BAM contain only unpaired reads?'...
       '<br><html>'];
  save_textfile(H,'report.html');
  return
end

% save all
fname = [individual '.dRanger_results.detail.all.mat'];
fprintf('Writing %s\n',fname);
save(fname,'X');
fname = [individual '.dRanger_results.detail.all.txt'];
fprintf('Writing %s\n',fname);
save_struct(X,fname);

% extract and save somatic rearrangements
idx = find(X.somatic_score>=MIN_FINAL_SOMATIC_SCORE);
X = reorder_struct(X,idx);
fname = [individual '.dRanger_results.detail.somatic.txt'];
fprintf('Writing %s\n',fname);
save_struct(X,fname);

% save simplified list of somatic rearrangements
flds = {'num','chr1','str1','pos1','chr2','str2','pos2','tumreads','normreads','class','span',...
  'site1','site2','fusion','quality','score','BPresult','BPsomaticratio','T_lenhomology','T_lenforeign','approxflag'};
Y = keep_fields(X,flds);
Y = orderfields(Y,flds);
Y.chr1 = chrlabel(Y.chr1);
Y.chr2 = chrlabel(Y.chr2);
strlabel = {'(+)';'(-)'};
Y.str1 = strlabel(Y.str1+1);
Y.str2 = strlabel(Y.str2+1);
Y.class = lower(Y.class);
fname = [individual '.dRanger_results.somatic.txt'];
fprintf('Writing %s\n',fname);
save_struct(Y,fname);

% write report.html
fprintf('Writing report.hmtl\n');
dRanger_write_report(individual,circospng,Y,MIN_FINAL_SOMATIC_SCORE);

