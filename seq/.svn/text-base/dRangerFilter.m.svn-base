function dRangerFilter(sample,P,X)
% dRangerFilter(sample,P[,X])
% can pass pre-loaded X to save time

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'results_name','dRanger_results');
P = impose_default_value(P,'perform_blast_filtering',true);

try

basedir = '/xchip/tcga_scratch/lawrence';
fname = [basedir '/' sample '/' P.results_name '.txt'];

if ~exist('X','var')
  if ~exist(fname,'file'), error('%s not found',fname); end
  X = load_struct(fname);
  X = make_numeric(X,{'chr1','chr2','min1','min2','max1','max2','str1','str2','tumreads','normreads'});
end

% SCATTEREDNESS FILTER

if 0, try
  fprintf('Annotating scatteredness\n');
  X = make_numeric(X,{'stdev1','stdev2'});
  d = [basedir '/' sample '/tumor_isz'];
  maxstdev = 200;
  setsizes = [2:15 20 30];
  p = dRanger_calculate_expected_stdev_distrib(d,maxstdev,setsizes);
  X.scatterP1 = nan(slength(X),1);
  X.scatterP2 = nan(slength(X),1);
  for r=1:slength(X)
    if ~mod(r,100), fprintf('%d/%d ',r,slength(X)); end
    for e=1:2
      if e==1, s = X.stdev1(r); else s = X.stdev2; end
      stdidx = max(1,min(maxstdev,round(s)));
      nr = X.tumreads(r);
      if nr<=setsizes(end), sszidx = find(setsizes<=nr,1,'last'); else sszidx = length(setsizes); end
      if e==1, X.scatterP1(r)=p(sszidx,stdidx); else X.scatterP2(r)=p(sszidx,stdidx); end
    end
  end
catch me
  fprintf('Scatteredness metric could not be added.  Did InsertSizeByLane complete?\n');
end, end

% ADD FILTERING METRICS FROM BOOM FILE
X.filterQ = zeros(slength(X),1);
X.filterW = zeros(slength(X),1);
t_boomdir = [basedir '/' sample '/tumor.boom'];
n_boomdir = [basedir '/' sample '/normal.boom'];
if exist(t_boomdir,'dir') && exist(n_boomdir,'dir')
%  try
    fprintf('Adding filtering metrics from boom file.\n');
    X = dRanger_add_filtering_from_boom_part2(X,t_boomdir,n_boomdir);
    X.filterQ = 1*(X.fmapqzN1>0.1 | X.fmapqzN2>0.1 | X.fmapqzT1>0.1 | X.fmapqzT2>0.1);
    lim = max(4,min(10,X.tumreads));
    X.filterW = 1*(X.nuwpN1>lim | X.nuwpN2>lim | X.nuwpT1>lim | X.nuwpT2>lim);
%  catch me
%    fprintf('    ...failed.\n');
%  end
else
  fprintf('Boom file not found: cannot retrieve filtering metrics.\n');
end

% BLAST FILTERING

if P.perform_blast_filtering
  fprintf('BLAST filtering...\n');
  X = dRanger_blast_filter(sample,P,X);
end

% COMBINE FILTERING CRITERIA

X = make_numeric(X,{'filterB','filterQ','filterW'});
X.filterBQW = 1*(X.filterB>0 | X.filterQ>0 | X.filterW>0);
X.filter = X.filterBQW;

% SORT AND SAVE

X = sort_struct(X,{'filter','normreads','tumreads'},[1 1 -1]);
fprintf('Saving results\n');
save_struct(X,fname);

fprintf('Done!\n');

catch me, excuse(me); end
