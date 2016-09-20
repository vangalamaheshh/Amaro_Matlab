function mutsig_preprocess(fhfile,isetname,outdir,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'build','*required');
P = impose_default_value(P,'build_dir','');
P = impose_default_value(P,'targlist','*required*');   %% list of targets to measure coverage on
P = impose_default_value(P,'categdir','*required*');   %% category set for measuring coverage (directory)
P = impose_default_value(P,'num_categories',4);        %% number of final categories to discover
P = impose_default_value(P,'mutation_blacklist',[]);
P = impose_default_value(P,'mutation_whitelist',[]);

P = impose_default_value(P,'consolidate_adjacent_muts_threshold',1);
P = impose_default_value(P,'bsub_queue','cga');

demand_file(fhfile);
demand_file(P.targlist);
demand_file(P.categdir);
if isfield(P,'categfile') && ~isempty(P.categfile)
  fprintf('P.categfile now ignored -- please remove\n');
end
P.categfile = [P.categdir '/all.fwb'];
demand_file(P.categfile);

Z = get_categs(P.categdir);
Z.num = str2double(Z.num);
if any(isnan(Z.num)) || any(Z.num<1), error('categs can''t have zero or nan'); end
mincateg = 1;
maxcateg = max(Z.num);

if isfield(P,'catlist') && ~isempty(P.catlist)
  error('Please remove P.catlist -- this will cause problems later with mutsig_run');
end

% check for possible build mismatches
if contains(P.categdir,'lawrence/db/.+')
  if (contains(P.build,'hg18') && contains(P.categdir,'hg19')) ||...
        (contains(P.build,'hg19') && (contains(P.categdir,'hg18') || ~contains(P.categdir,'hg19')))
    error('possible mismatch between P.build=%s and P.categdir=%s',P.build,P.categdir);
end,end
if (contains(P.build,'hg18') && contains(P.targlist,'hg19')) ||...
      (contains(P.build,'hg19') && contains(P.targlist,'hg18'))
  error('possible mismatch between P.build=%s and P.targlist=%s',P.build,P.targlist);
end

if isempty(P.build_dir)
  fprintf('Using hard-coded path to UCSC build dirs\n');
  P.build_dir = ['/xchip/cga1/annotation/db/ucsc/' P.build];
end
if ~exist(P.build_dir,'dir')
  fprintf('WARNING: Could not find %s\n', P.build_dir);
end

try
  ReferenceInfoObj.init(P.build_dir);
catch me
  fprintf('WARNING: failed to initialize ReferenceInfoObj\n');
end

% CREATE OUTPUT DIRECTORY
ensure_dir_exists(outdir);
outstem = [outdir '/' isetname];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PATIENT LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Preprocessing patient list\n');

pat = load_struct(fhfile);
flds = {'Individual_Id','maf_file_capture','somatic_mutation_coverage_capture','indel_maf_file_capture'};
demand_fields(pat,flds);
pat.name = pat.Individual_Id;
pat.maf = pat.maf_file_capture;
pat.wig = pat.somatic_mutation_coverage_capture;
pat.indel_maf = pat.indel_maf_file_capture;
pat.dataset = repmat({'capture'},slength(pat),1);

for i=1:slength(pat)
  if iscell(pat.wig{i}), error('Multiple WIGs per patient not yet implemented in this function'); end
  if iscell(pat.maf{i}), error('Multiple MAFs per patient not yet implemented in this function'); end
  if iscell(pat.indel_maf{i}), error('Multiple indels MAFs per patient not yet implemented in this function'); end
end

demand_file(pat.wig);
demand_file(pat.maf);
demand_file(pat.indel_maf);

patient_list = [outstem '.patients.txt'];
save_struct(pat,patient_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MUTATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Preprocessing mutations\n');

x = [];
for i=1:slength(pat)
  fprintf('Loading %s\n',pat.maf{i});
  x{end+1,1} = load_struct(pat.maf{i});
  x{end}.patient = repmat({pat.name{i}},slength(x{end}),1);
  x{end}.dataset = repmat({pat.dataset{i}},slength(x{end}),1);
  fprintf('Loading %s\n',pat.indel_maf{i});
  x{end+1,1} = load_struct(pat.indel_maf{i});
  x{end}.patient = repmat({pat.name{i}},slength(x{end}),1);
  x{end}.dataset = repmat({pat.dataset{i}},slength(x{end}),1);
end
x = concat_structs_keep_all_fields(x);

x = preprocess_mutations(x,P.categdir,P);

% save final maf
mutation_list = [outstem '.maf'];
save_struct(x,mutation_list);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  COVERAGE AND CATEGORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Preprocessing coverage\n');

% extract coverage from wigs
covdir = [outstem '.covfiles']; ensure_dir_exists(covdir);
pat.covfile = regexprep(pat.name,'(.*)',[covdir '/$1.somatic_coverage.txt']);
PP=[]; PP.bsub_queue = P.bsub_queue;
extract_from_wig(pat.name,P.targlist,pat.covfile,P.categfile,P.build,P.build_dir,mincateg,maxcateg,pat.wig,PP);
% (submits jobs to LSF and waits for them to finish)

% gather
coverage_gather(patient_list,covdir,outdir,isetname,P.targlist,P.categdir,P.build,P.build_dir,mutation_list,P.num_categories,P)

fprintf('MutSig preprocessing completed successfully!\n');
