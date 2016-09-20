function fh_dRangerRun(individual,tdir,ndir,tminmapq,minpairs,windowsize,nminwindow,normpaneldb,minsomratio,nminspanfrac,minscoreforbp,min_ignore_matchingnormal,build)
% fh_dRangerRun(individual,tdir,ndir,tminmapq,minpairs,windowsize,nminwindow,normpaneldb,minsomratio,nminspanfrac,minscoreforbp,build)
%
% Mike Lawrence 2010

if nargin<3, error(['Usage: fh_dRangerRun(individual,tdir,ndir,tminmapq,minpairs,windowsize,'...
    'nminwindow,normpaneldb,minsomratio,nminspanfrac,minscoreforbp,min_ignore_matchingnormal,build)']); end

if ~exist('tminmapq','var'), tminmapq = 5; end
if ~exist('minpairs','var'), minpairs = 2; end
if ~exist('windowsize','var'), windowsize = 2000; end
if ~exist('nminwindow','var'), nminwindow = 2000; end
if ~exist('normpaneldb','var'), normpaneldb = ''; end
if ~exist('minsomratio','var'), minsomratio = 50; end
if ~exist('nminspanfrac','var'), nminspanfrac = 0.5; end
if ~exist('minscoreforbp','var'), minscoreforbp = 0.01; end
if ~exist('min_ignore_matchingnormal','var'), min_ignore_matchingnormal = Inf; end
if ~exist('build','var'), build = 'hg18'; end

% note: build can be either hg18/hg19/etc.  OR  the absolute path to a Refseq matfile for load_refseq

if ~isnumeric(tminmapq), tminmapq = str2double(tminmapq); end
if ~isnumeric(minpairs), minpairs = str2double(minpairs); end
if ~isnumeric(windowsize), windowsize = str2double(windowsize); end
if ~isnumeric(nminwindow), nminwindow = str2double(nminwindow); end
if ~isnumeric(minsomratio), minsomratio = str2double(minsomratio); end
if ~isnumeric(nminspanfrac), nminspanfrac = str2double(nminspanfrac); end
if ~isnumeric(minscoreforbp), minscoreforbp = str2double(minscoreforbp); end
if ~isnumeric(min_ignore_matchingnormal), min_ignore_matchingnormal = str2double(min_ignore_matchingnormal); end

fprintf('fh_dRangerRun\n');
fprintf('  individual = %s\n', individual);
fprintf('  tdir = %s\n', tdir);
fprintf('  ndir = %s\n', ndir);
fprintf('  tminmapq = %d\n', tminmapq);
fprintf('  minpairs = %d\n', minpairs);
fprintf('  windowsize = %d\n', windowsize);
fprintf('  nminwindow = %d\n', nminwindow);
fprintf('  normpaneldb = %s\n', normpaneldb);
fprintf('  minsomratio = %d\n', minsomratio);
fprintf('  nminspanfrac = %0.2f\n', nminspanfrac);
fprintf('  minscoreforbp = %0.2f\n', minscoreforbp);
fprintf('  min_ignore_matchingnormal = %0.2f\n', min_ignore_matchingnormal);
fprintf('  build = %s\n', build);

P=[];
P = impose_default_value(P,'dRanger_input_matfile','all.weird.pairs.mat');
P = impose_default_value(P,'dRanger_input_iszfile','all.isz');
P = impose_default_value(P,'weirdpair_log10_threshold',-3.5);
P = impose_default_value(P,'tumor_mapping_quality_cutoff',tminmapq);
P = impose_default_value(P,'minpairs',minpairs);
P = impose_default_value(P,'discovery_window_size',windowsize);
P = impose_default_value(P,'minimum_normal_window',nminwindow);
P = impose_default_value(P,'minimum_normal_span_frac',nminspanfrac);
P = impose_default_value(P,'panel_of_normals_database',normpaneldb);
P = impose_default_value(P,'minsomratio',minsomratio);
P = impose_default_value(P,'results_name','dRanger_results');
P = impose_default_value(P,'min_score_to_perform_BreakPointer',minscoreforbp);
P = impose_default_value(P,'min_ignore_matchingnormal',min_ignore_matchingnormal);
P = impose_default_value(P,'build',build);
P = impose_default_value(P,'save_TT_and_NN_in_output_matfile',false);
P = impose_default_value(P,'delete_input_matfiles_on_success',true);

disp(P)

% LOAD INPUT DATA
fprintf('Loading input data\n');

% load insert-size distributions
iszT = load_isz_file([tdir '/' P.dRanger_input_iszfile]);
iszN = load_isz_file([ndir '/' P.dRanger_input_iszfile]);

% load weird-pair data
[T TT] = load_dRanger_input([tdir '/' P.dRanger_input_matfile],P);
[N NN] = load_dRanger_input([ndir '/' P.dRanger_input_matfile],P);

% check for case of "no data"
if ~isempty(T) && ~isempty(N) && iszT.tot>0 && iszN.tot>0

  % if necessary, increase stringency of cutoff for tumor weird-pair status
  % (depending on tumor+normal insert-size distributions)
  t_wp_cutoff = dRanger_calculate_weirdpair_cutoff_from_isz(iszT, P.weirdpair_log10_threshold);
  n_wp_cutoff = dRanger_calculate_weirdpair_cutoff_from_isz(iszN, P.weirdpair_log10_threshold);
  wp_cutoff = max(t_wp_cutoff,n_wp_cutoff);
  tspan = T(:,8)-T(:,4); tspan(T(:,2)~=T(:,6)) = nan;
  idx = find(tspan<wp_cutoff);
  if ~isempty(idx)
    fprintf('Removing %d/%d tumor readpairs that have spans within the 10^%.1f cutoff\n',...
      length(idx),size(T,1),P.weirdpair_log10_threshold);
    TT = reorder_struct_exclude(TT,idx);
    T(idx,:) = [];
  end
%  n_cutoff = dRanger_calculate_weirdpair_cutoff_from_isz(iszN, P.normal_weirdpair_log10_threshold);
%  nspan = N(:,8)-N(:,4); nspan(N(:,2)~=N(:,6)) = nan;
%  idx = find(nspan<n_cutoff);
%  if ~isempty(idx)
%    fprintf('Removing %d/%d normal reads that have spans within the 10^%.1f cutoff\n',...
%      length(idx),size(N,1),P.weirdpair_log10_threshold);
%    NN = reorder_struct_exclude(NN,idx);
%    N(idx,:) = [];
%  end

  % filter tumor data to remove low-quality reads
  idx = find(TT.qual1>=P.tumor_mapping_quality_cutoff & TT.qual2>=P.tumor_mapping_quality_cutoff);
  TT = reorder_struct(TT,idx);
  T = T(idx,:);
  
  % RUN CORE ALGORITHM
  [X TT.ridx NN.ridx] = dRanger_core_algorithm(T,N,P);
  
  % ADD CLASS
  X.class = classify_rearrangements(X);

  % ORDER FIELDS
  X.individual = repmat({individual},slength(X),1);
  flds = {'individual','num','chr1','str1','pos1','chr2','str2','pos2','class','span','tumreads','normreads'};
  X = orderfields_first(X,flds);

  % SCREEN AGAINST PANEL OF NORMALS (if specified)
  if ~isempty(P.panel_of_normals_database)
    X = dRanger_screen_against_panel_of_normals(X,P);
    flds = [flds 'normpanelreads','normpanelsamps','normpaneldetails'];
    X = orderfields_first(X,flds);
  end
  
  % ANNOTATE SITES BY GENE AND FUSION STATUS
  X = dRanger_annotate_sites(X,P);
  
  % SCORING
  
  % extract nuwp and fmapqz
  fprintf('Adding filtering metrics from *.fmapqz and *.nuwp files.\n');
  X = dRanger_add_filtering_from_boom_part2(X,tdir,ndir,struct('add_coverage',false,'add_avgnmm',false));
  
  % evaluate standard deviation of starting position of supporting reads
  std_distrib = dRanger_calculate_expected_stdev_distrib(iszT);
  mi = nan(slength(X),1);
  for i=1:slength(X)
    mi(i) = find(std_distrib(:,1)<=X.tumreads(i),1,'last');
  end
  
  X.zstdev1 = (X.stdev1 - std_distrib(mi,2)) ./ std_distrib(mi,3);
  X.zstdev2 = (X.stdev2 - std_distrib(mi,2)) ./ std_distrib(mi,3);
  
  % (range z-score no longer used; provides little new information beyond stdev z-score)
  % X.zrange1 = (X.range1 - std_distrib(mi,4)) ./ std_distrib(mi,5);
  % X.zrange2 = (X.range2 - std_distrib(mi,4)) ./ std_distrib(mi,5);

  % compute quality
  X.quality = dRanger_calculate_quality(X);

  % compute score
  X.score = X.tumreads .* X.quality;
  X.score(round(X.score)<P.minpairs) = 0;     % impose threshold = minpairs

  % compute somatic score
  if isfield(X,'normpanelreads')
      normreads = max(X.normreads,X.normpanelreads);
      ignore_matchingnormal = (X.tumreads ./ X.normreads) >= P.min_ignore_matchingnormal;
      if any(ignore_matchingnormal)
        normreads(ignore_matchingnormal) = X.normpanelreads(ignore_matchingnormal);
      end
  else
      normreads = X.normreads;
  end
  X.somatic = (X.tumreads ./ normreads) >= P.minsomratio;
  X.somatic_score = X.score .* X.somatic;

  % sort rearrangements by somatic score 
  a = nan(slength(X),1);
  [X ord] = sort_struct(X,'somatic_score',-1);
  a(ord) = 1:slength(X);
  a = as_column(a);
  TT.ridx = nansub(a,TT.ridx);
  NN.ridx = nansub(a,NN.ridx);
    
  % coarse-filter results for BreakPointer
  X.BPtry = (X.score>=P.min_score_to_perform_BreakPointer);
  Xbp = reorder_struct(X,X.BPtry);

else    % NO DATA
  X = []; X.NO_DATA = [];
  Xbp = X;
  TT = [];
  NN = [];
end


% SAVE list of rearrangements for BreakPointer to assemble
fname = [individual '.' P.results_name '.forBP.txt'];
fprintf('Saving %s\n',fname);
save_struct(Xbp,fname);

% SAVE results to matfile
fname = [individual '.' P.results_name '.mat'];
fprintf('Saving %s\n',fname);
if P.save_TT_and_NN_in_output_matfile
    save(fname,'X','TT','NN','-v7.3');
else
    save(fname,'X','-v7.3');
end

% DELETE input matfiles
if P.delete_input_matfiles_on_success
  try
    system(['rm ' tdir '/' P.dRanger_input_matfile]);
    system(['rm ' ndir '/' P.dRanger_input_matfile]);
  catch me
    fprintf('Error during attempt to delete input matfiles at end of fh_dRangerRun\n');
    disp me
    disp me.message
  end
end
