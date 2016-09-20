function dRangerRun(sample,P)
% dRangerRun(sample,P)
%

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR1');
P = impose_default_value(P,'dRanger_input_matfile','dRanger_input.mat');
P = impose_default_value(P,'flip_TN',false);
P = impose_default_value(P,'remove_duplicates',true);
P = impose_default_value(P,'tumor_mapping_quality_cutoff',5);
P = impose_default_value(P,'min_pairs',2);
P = impose_default_value(P,'discovery_window_size',2000);
P = impose_default_value(P,'minimum_normal_window',400);
P = impose_default_value(P,'pause_after_core_algorithm',false);
P = impose_default_value(P,'save_X_TT_and_NN',true);
P = impose_default_value(P,'results_name','dRanger_results');
P = impose_default_value(P,'add_annotation',true);

if isfield(P,'window_size')
  fprintf('Note: "window_size" parameter has been changed to "discovery_window_size"\n');
  P.discovery_window_size = P.window_size;
end

try

tic

fprintf('dRangerRun\n  sample = %s\n',sample);

basedir = '/xchip/tcga_scratch/lawrence';
fname = [basedir '/' sample '/' P.results_name '.txt'];
if exist(fname,'file')
  fprintf('Results file already exists.  Delete it if you want to re-run dRanger.\n');
  fprintf('Type "return" to continue on to check that filtering is complete.');
  keyboard
  fprintf('Checking to see that filtering is complete.\n');
  dRangerFilter(sample,P);
  return
end  

% LOAD INPUT DATA
fprintf('Loading input data\n');
fsuffix = [P.dRangerPreprocess_output_dir_suffix '/' P.dRanger_input_matfile];
% load and de-dup
[T TT] = load_dRanger_input([basedir '/' sample '/tumor_' fsuffix]);
[N NN] = load_dRanger_input([basedir '/' sample '/normal_' fsuffix]);
% flip T/N if requested
if P.flip_TN, tmp = TT; TT = NN; NN = tmp; tmp = T; T = N; N = tmp; end
% filter tumor data to remove low-quality reads
idx = find(TT.qual1>=P.tumor_mapping_quality_cutoff & TT.qual2>=P.tumor_mapping_quality_cutoff);
TT = reorder_struct(TT,idx);
T = T(idx,:);

% RUN CORE ALGORITHM
%[trans,F,E,filtered_E] = dRanger_core_algorithm(T,N,P);
%X = dRanger_reformat(trans);
[X TT.ridx NN.ridx] = dRanger_core_algorithm(T,N,P);
if P.pause_after_core_algorithm
  fprintf('Pausing after core algorithm\n');
  keyboard
end

% ADD CLASS
X.class = classify_rearrangements(X);

% ANNOTATE GENES
if P.add_annotation
  X = dRanger_annotate(X);                       % add gene/site annotations
end

% ANNOTATE VALDATION (not for all samples)
if strcmp(sample,'gbm/0188/wgs')
  X = dRanger_xreftoval(X);                      % cross-reference to GBM0188 validation data
end

% SAVE
if P.save_X_TT_and_NN  % default = true
  fname = [basedir '/' sample '/' P.results_name '_X_TT_NN.mat'];
  fprintf('Saving %s\n',fname);
  save(fname,'X','TT','NN','-v7.3');   % ~3 min
end

fname = [basedir '/' sample '/' P.results_name '.txt'];
fprintf('Saving %s\n',fname);
save_struct(X,fname);

% PROCEED TO FILTERING STEP
dRangerFilter(sample,P,X);

% DONE
toc

catch me, excuse(me); end
