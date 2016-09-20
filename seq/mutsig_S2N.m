function mutsig_S2N(P)

if nargin<1, error('need P struct'); end
P = impose_default_value(P,'isetname','dataset');
P = impose_default_value(P,'build','*required*');
P = impose_default_value(P,'mutfile','*required*');
P = impose_default_value(P,'covfile','*required*');
P = impose_default_value(P,'outdir','*required*');
P = impose_default_value(P,'outfile',[]);
P = impose_default_value(P,'psubsets',[]);         % default will be:  one subset consisting of all patients

% load data
M = mutsig_load(P);

% error-check patient subsets
if isempty(P.psubsets), P.psubsets = {M.patient.name}; end
if ~iscell(P.psubsets), error('P.psubsets should be a cell array'); end
for i=1:length(P.psubsets)
  p = P.psubsets{i};
  if ~iscellstr(p), error('P.psubsets should be a cell array of cellstrs'); end
  z = listmap(p,M.patient.name);
  if any(isnan(z))
    disp(M.patient.name(isnan(z)));
    fprintf('WARNING: no such patient(s)\n');
    P.psubsets{i} = P.psubsets{i}(~isnan(z));
  end
end
if ~isempty(P.outfile)
  if ~iscellstr(P.outfile), error('P.outfile should be a cellstr array'); end
  if length(P.outfile)~=length(P.psubsets), error('P.outfile should be same length as P.psubsets'); end
  if length(unique(P.outfile))~=length(P.psubsets), error('P.outfile should be list of unique filenames'); end
  if ~isempty(grep('^/',P.outfile)), error('P.outfile should not be absolute paths--use P.outdir'); end
else
  P.outfile = cell(length(P.psubsets),1);
  for i=1:length(P.psubsets)
    P.outfile{i} = ['sig_genes.' num2str(i) '.txt'];
  end
end 

% find bagels (in principled way) if not using precomputed ones
if ~isfield(M.gene,'bagel')
  P = impose_default_value(P,'min_neighbors',5);
  P = impose_default_value(P,'max_neighbors',30);
  M = neighborhood_mutsig_fromM_v2(M,P);
end

% open write directory
ede(P.outdir);

% load COSMIC
COS = load_cosmic_database(P);

% for each sample subset, run the significance calculation and save results
for i = 1:length(P.psubsets)
  fprintf('Subset %d/%d\n',i,length(P.psubsets));
  P.signal_pat_subset = P.psubsets{i};
  M = simple_calc_rates_v3(M,P);
  stem = [P.outdir '/' P.outfile{i}];
  fname = [stem '.sig_genes.txt'];
  save_struct(M.g,fname);
  fname = [stem '.sig_genes.mat'];
  z = M.gene; save(fname,'z');
  P.cosmic_report_filename = [stem '.cosmic.sig_genes.txt'];
  P.cosmic_report2_filename = [stem '.cosmic.sig_genes2.txt'];
  P.cosmic_mutations_outname = [stem '.cosmic.mutations.txt'];
  M.mutrate.tot.hat = M.rnon_c(end);
  P.patients_to_include = P.psubsets{i};
  perform_COSMIC_overlap_analysis(M,P,COS)
end




