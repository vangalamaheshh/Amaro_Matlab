function get_global_coverage_stats(samples,outname,filter,P)
% get_coverage_global_stats(samples,outname,filter)
%
% samples
% "CLASSICAL MODE"
%    e.g.    gbm/0188/wgs
%            cll/je/wgs
%            ov/0725/c2k
%            ov/0751/c6k
%    (can be a cell array of strings)
% "FULL-PATH MODE"
%    samples should be a two-column cell array,
%          where column 1 lists the CoverageByBase output diretory for the tumor
%          and column 2 is the normal
%    outname is a one-column cell array listing the full path of the output file
%          (a temporary directory will also be created here)
%
% outname
%    e.g.    coverage_by_zone.txt
%
% filter
%    e.g.    context
%            zone
%            c2k 
%            c6k
%    (multiple filters not currently supported)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'covbb_file_extension','cbb');
P=impose_default_value(P,'mincateg',[]);
P=impose_default_value(P,'maxcateg',[]);

F = load_struct(['/xchip/cga1/lawrence/db/' filter '/categs.txt'],'%f%s');
if ~issorted(F.num), error('F.num should be sorted'); end

if isempty(P.mincateg) || isempty(P.maxcateg)
  % infer from categs.txt
  P.mincateg = F.num(1);
  P.maxcateg = F.num(end);
end

minmax = [num2str(P.mincateg) ' ' num2str(P.maxcateg)];
ncateg = P.maxcateg-P.mincateg+1;

if ~exist('filter','var'), filter=[]; end
if isempty(filter), error('Non-filtered operation not currently supported'); end
if ~iscell(samples), samples = {samples}; end
nsamp = size(samples,1);

if F.num(1)~=P.mincateg, error('F.num(1) doesn''t match P.mincateg'); end
if F.num(end)~=P.maxcateg, error('F.num(end) doesn''t match P.maxcateg'); end

scriptname = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/get_coverage_stats.pl';
basedir = '/xchip/cga1/lawrence';
if ~isempty(filter), filterdir = ['/xchip/cga1/lawrence/db/' filter];
else filterdir = 'none'; end

% establish outfile{} tumdir{} normdir{}
tumdir = cell(nsamp,1);
normdir = cell(nsamp,1);
outfile = cell(nsamp,1);
if size(samples,2)==1       % "CLASSICAL-MODE"
  for i = 1:nsamp
    sample = samples{i};
    direc = ['/xchip/cga1/lawrence/' sample];
    tumdir{i} = [direc '/tumor_cbb'];
    normdir{i} = [direc '/normal_cbb'];
    outfile{i} = [direc '/' outname];
  end
elseif size(samples,2)==2   % "FULL-PATH MODE"
  tumdir = samples(:,1);
  normdir = samples(:,2);
  outfile = outname;
else
  error('Unknown input format');
end

% first see if final output file exists for any of these samples
already_done = false(nsamp,1);
for i=1:nsamp
  dout = dir(outfile{i});
  uptodate = false;
  if ~isempty(dout)
    uptodate = true;
    for c=1:24
      tumfile = [tumdir{i} '/chr' num2str(c) '.' P.covbb_file_extension];
      normfile = [normdir{i} '/chr' num2str(c) '.' P.covbb_file_extension];
      dtum = dir(tumfile);
      dnorm = dir(normfile);
      if ~isempty(dtum) && dtum.datenum>dout.datenum, uptodate = false; end;
      if ~isempty(dnorm) && dnorm.datenum>dout.datenum, uptodate = false; end;
    end
    if uptodate, already_done(i) = true; end
  end
end
samples = samples(~already_done,:);
outfile = outfile(~already_done);
nsamp = size(samples,1);
if nsamp==0, fprintf('All files are up-to-date\n'); return; end

% make sure all needed files exist

all_ok = true;
for i=1:nsamp
    for c=1:24
      tumfile = [tumdir{i} '/chr' num2str(c) '.' P.covbb_file_extension];
      normfile = [normdir{i} '/chr' num2str(c) '.' P.covbb_file_extension];
      d = dir(tumfile);
      if isempty(d) | isempty(d.bytes)
        fprintf('%s: file not found\n',tumfile);
        all_ok = false;
      end
      d = dir(normfile);
      if isempty(d) | isempty(d.bytes)
        fprintf('%s: file not found\n',normfile);
        all_ok = false;
end,end,end

if ~strcmp('none',filterdir)
  for c=1:24
  filterfile = [filterdir '/chr' num2str(c) '.txt'];
  d = dir(filterfile);
  if isempty(d) | isempty(d.bytes)
    fprintf('%s: file not found\n',filterfile);
    all_ok = false;
end,end,end

if ~all_ok, error('Some required files not found.  Quitting.\n'); end

% send out all jobs

all_done = false;
while(~all_done)
  jobs = [];
  cmds = {}; banners = {};
  all_done = true;
  for i=1:nsamp
    outdir = [outfile{i} '_temp'];
    if ~exist(outdir,'dir'), mkdir(outdir); end
    for c=1:24
      tumfile = [tumdir{i} '/chr' num2str(c) '.' P.covbb_file_extension];
      normfile = [normdir{i} '/chr' num2str(c) '.' P.covbb_file_extension];
      dtum = dir(tumfile);
      dnorm = dir(normfile);
      statfile = [outdir '/chr' num2str(c) '.stats'];
      if exist(statfile,'file')
        dstat = dir(statfile);
        if ~isempty(dstat) && ~isempty(dstat.bytes)
          if dstat.datenum>=dtum.datenum && dstat.datenum>=dnorm.datenum, continue;
          else fprintf('%s needs to be refreshed:\n',statfile); end
        end
      end
      all_done = false;
      banners{end+1,1} = [num2str(i) 'GLOB' num2str(c)];
      cmds{end+1,1} = ['"perl ' scriptname ' ' tumdir{i} ' ' normdir{i} ' ' filterdir ' ' outdir ' ' num2str(c) ' ' minmax '"'];
      % submit jobs as they accumulate
      if length(cmds)>=60
        jobs = [jobs; bsub(cmds,banners)];
        cmds = {}; banners = {};
      end
    end
  end
  if ~isempty(cmds)
    jobs = [jobs; bsub(cmds,banners)];
    cmds = {}; banners = {};
  end
  if (~all_done) bwait(jobs); end
end

% integrate statistics across chromosomes
fprintf('Integrating statistics across chromosomes:\n');

for i=1:nsamp, fprintf('Sample %d/%d\n',i,nsamp);
  try
    rlen = getBAMFileReadLength(samples{i});
  catch me
    rlen = 101;
  end
  outdir = [outfile{i} '_temp'];
  X=F;
  z = zeros(ncateg,1); X.terr=z; X.tseqbp=z; X.nseqbp=z; X.tseqreads=z; X.nseqreads=z;
  X.tdepth=z; X.ndepth=z; X.callablebp=z; X.callablefrac=z;
  for c=1:24
    statfile = [outdir '/chr' num2str(c) '.stats'];
    tmp = load_struct(statfile,'%f%f%f%f%f%f',0);
    tmp = rename_fields(tmp,{'col1','col2','col3','col4','col5','col6'},...
    {'chr','categ','terr','tseqbp','nseqbp','callablebp'});
    if any(tmp.categ~=X.num), error('tmp.categ doesn''t match X.num'); end
    X.terr=X.terr+tmp.terr;
    X.tseqbp=X.tseqbp+tmp.tseqbp;
    X.nseqbp=X.nseqbp+tmp.nseqbp;
    X.callablebp=X.callablebp+tmp.callablebp;
  end
  % add row of summary statistics
  X.num(end+1) = nan;
  X.name{end+1} = 'Total';
  X.terr(end+1) = sum(X.terr);
  X.tseqbp(end+1) = sum(X.tseqbp);
  X.nseqbp(end+1) = sum(X.nseqbp);
  X.callablebp(end+1) = sum(X.callablebp);
  % compute derived columns
  X.tseqreads = round(X.tseqbp/rlen);
  X.nseqreads = round(X.nseqbp/rlen);
  X.tdepth = X.tseqbp ./ X.terr;
  X.ndepth = X.nseqbp ./ X.terr;
  X.callablefrac = X.callablebp ./ X.terr;
  % write output file
  save_struct(X,outfile{i});
end


