function dRangerRun(sample,P)
% dRangerRun(sample,P)
%
% old parameter style: (sample,min_pairs,window_size,C_threshold,H_threshold,L_threshold)

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR');
P = impose_default_value(P,'dRangerPrepare_output_file','all.weird.joined.mat');
P = impose_default_value(P,'results_name',[]);

% eventually will be:
%P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
%P = impose_default_value(P,'dRangerPrepare_output_file','stringent_pairs.mat');

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'min_pairs',4);
P = impose_default_value(P,'window_size',2000);
P = impose_default_value(P,'C_threshold',2000);
P = impose_default_value(P,'H_threshold',0.9);
P = impose_default_value(P,'L_threshold',1.1);
P = impose_default_value(P,'Z_threshold',0.3);
P = impose_default_value(P,'save_matfile',false);
P = impose_default_value(P,'pause_after_core_algorithm',false);

try

tic

direc = ['/xchip/tcga_scratch/lawrence/' sample];

[T N] = dRangerLoad(sample,P);
[trans,F,E,filtered_E] = dRanger2(T,N,P.min_pairs,P.window_size,[]);
X = dRanger_reformat(trans);
X = sort_struct(X,{'normreads','tumreads'},[1 -1]);

if P.pause_after_core_algorithm
  fprintf('Pausing after core algorithm\n');
  keyboard
end

% GATHER OTHER DATA
X = dRanger_annotate(X);                       % add gene/site annotations
X = dRanger_annotate_centel(X,1000000,25000);  % annotate centromeres,telomeres
X = dRanger_xreftoval(X);                      % cross-reference to GBM0188 validation data

% FILTER RESULTS
fprintf('Filtering results...\n');
if P.cancer_sample, ts = 'tumor'; else ts = 'sample'; end
[X.count1 X.count2] = dRanger_find_local_coverage(X,[direc '/' ts '_' P.dRangerPreprocess_output_dir_suffix '/']);
X = dRanger_find_densities(X,T);
X.H = 2*X.tumreads./(X.tot1tum+X.tot2tum);
X.L = X.Lsupp1tum./X.tumreads;
X.filterC = 1*(X.count1>P.C_threshold | X.count2>P.C_threshold);
X.filterH = 1*(X.H<P.H_threshold);
X.filterL = 1*(X.L>P.L_threshold);
% X.filterHCL = 1*(X.filterH|X.filterC|X.filterL);
[X.mqz1 X.mqz2] = dRanger_find_local_MQZ_counts(X,[direc '/' ts '_' P.dRangerPreprocess_output_dir_suffix '/']);
X.Z1 = X.mqz1 ./ X.count1;
X.Z2 = X.mqz2 ./ X.count2;
X.filterZ = 1*(X.Z1>P.Z_threshold|X.Z2>P.Z_threshold);

% SAVE RESULTS
fprintf('Saving results...\n');
name2 = upper(regexprep(sample,'/','-'));
if isempty(P.results_name), results_name = [name2 '_dRanger_results'];
else results_name = P.results_name; end
save_struct(X,[direc '/' results_name '.txt']);
if P.save_matfile
  fprintf('Saving matfile...\n');
  matfile = [direc '/' results_name '.mat'];
  save(matfile,'-v7.3');
end

% BLAST FILTERING
dRanger_blast_filter(sample,P);

% DONE
toc

catch me, excuse(me); end
