function coverage_gather(patient_list,covdir,outdir,individual_set_id,target_list,categdir,build,build_dir,mutation_list,num_categories,P)
% called by mutsig_preprocess.m (outside Firehose)
%  or fh_MutSigCoverageGather.m (inside Firehose)

if nargin<10, error('not enough arguments'); end

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'add_double_null_category',true);
P=impose_default_value(P,'add_indel_null_category',true);
P=impose_default_value(P,'add_indel_category',false);
P=impose_default_value(P,'add_null_category',false);
P=impose_default_value(P,'impute_full_coverage',false);

if P.add_indel_null_category && (P.add_indel_category || P.add_null_category)
  error('indel+null category is incompatible with indel, null categories');
end

outstem = [outdir '/' individual_set_id];

categ_fwb = fullfile(categdir, 'all.fwb');
categ_list = fullfile(categdir, 'categs.txt');
demand_file({categ_fwb;categ_list});

if ischar(patient_list)
  pat = load_struct(patient_list);
elseif isstruct(patient_list)
  pat = patient_list;
  demand_fields(pat,'name');
else
  error('invalid format for patient_list parameter');
end

if ~P.impute_full_coverage
  % make sure coverage files exist
  if ~isfield(pat,'fwbfile')
    pat.fwbfile = regexprep(pat.name,'(.*)',[covdir '/$1.somatic_coverage.fwb']);
  end
  demand_files(pat.fwbfile);
  if ~isfield(pat,'covfile')
    pat.covfile = regexprep(pat.name,'(.*)',[covdir '/$1.somatic_coverage.txt']);
  end
  tmp = regexprep(pat.covfile,'(.*)','$1.bin');
  binmode = true;
  for i=1:length(tmp), if ~exist(tmp{i},'file'), binmode = false; break; end, end
  if binmode, pat.covfile = tmp; end
  demand_files(pat.covfile);
  
  % create summed FWB file
  summed_fwb = [outstem '.summed_coverage.fwb'];
  create_summed_FWB_coverage_file(pat.fwbfile,target_list,summed_fwb,P);
end

% categories
flds = {'left','from','change','right','autoname','name','type'};

if ischar(num_categories)
  if num_categories(1)=='/'
    % user-specified category set
    categfile = num_categories;
    fprintf('Skipping automatic category discover.\nInstead, using categories file "%s".\n',categfile);
    demand_file(categfile);
    K = load_struct(categfile);
    demand_fields(K,flds);
  else
    num_categories = str2double(num_categories);
  end
end

if isnumeric(num_categories)
  % automatic category discovery
  if num_categories<1, error('num_categories must be at least 1'); end
  if num_categories>6
    fprintf('WARNING: num_categories>6:  may take a very long time to compute best category set\n');
  end
  if ~P.impute_full_coverage
    N = get_covered_territory_tally_from_FWB(summed_fwb,categ_fwb);
  else
    %% get category tally of territory (this is a hacky fix, because it duplicates what happens in load_coverage_compactly...)
    %% TO DO: start building "C" object earlier (i.e. here), and do the territory calculation up front.
    T = load_target_file(target_list);
    nt = slength(T);
    cat = load_struct(categ_list);
    if strcmp(cat.num{1},'0'), error('Category list includes a zero category!'); end
    ncat = slength(cat);
    fprintf('Tabulating territory breakdown by categories.\n');
    N = zeros(ncat,1);
    fwb = org.broadinstitute.cga.tools.seq.FixedWidthBinary(categ_fwb);
    for ti=1:nt, if ~mod(ti,10000), fprintf('%d/%d ',ti,nt); end
      k = fwb.get(T.chr(ti),T.start(ti),T.end(ti));
      N = N + histc(k,1:ncat);
    end, fprintf('\n');
    fwb.close();
    N = N * slength(pat);
  end
  N = collapse_categories_to_65(N,categ_list);   % in case context is a more complicated one, e.g. c65e
  if ~any(N>0), error('Something went wrong with coverage processing: no usable coverage data!'); end
  m = load_struct(mutation_list);
  m = reorder_struct(m,strcmp('SNP',m.classification));
  bases = {'A';'C';'G';'T'};
  m.newbase = find_newbase(m);
  m.to = listmap(m.newbase,bases);
  n = hist2d_fast(str2double(m.context65),m.to,1,65,1,4);
  Nn = collapse_Nn_65_to_32([N n]);
  PP=[]; PP.max_k = num_categories;
  PP.mutcategs_report_filename = [outstem '.mutcateg_discovery.txt'];
  Ks = find_mut_categs(Nn,PP);
  K = Ks{num_categories};
  K = keep_fields(K,flds);
  if P.add_indel_category
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'indel'};
    K2.name = {'indel'};
    K2.type = {'non-point'};
    K = concat_structs({K,K2});
  end
  if P.add_null_category
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'null'};
    K2.name = {'null'};
    K2.type = {'non-point'};
    K = concat_structs({K,K2});
  end
  if P.add_indel_null_category 
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'indel+null'};
    K2.name = {'indel+null'};
    K2.type = {'non-point'};
    K = concat_structs({K,K2});
  end
  if P.add_double_null_category
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'double_null'};
    K2.name = {'double_null'};
    K2.type = {'non-point'};
    K = concat_structs({K,K2});
  end
end

% save categories
fname = [outstem '.mutcategs.txt'];
fprintf('Saving %s\n',fname);
save_struct(K,fname);

% process coverage
C = []; C.sample = pat; C.ns = slength(C.sample);
C.file.targ = target_list;
C.file.categdir = categdir;
C1 = load_coverage_compactly_and_convert_to_mat(C,K,P);
fname = [outstem '.coverage.mat'];
fprintf('Saving %s\n',fname);
save(fname,'C1','-v7.3');

