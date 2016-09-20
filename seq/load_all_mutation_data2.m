function M = load_all_mutation_data2(P,isetname,build)
% Mike Lawrence 2010-01-27

if ischar(P)
  indir = P;
  P=[];
  P.mutfile = decell(direc([indir '/*.final_analysis_set.maf']));
  P.covfile = decell(direc([indir '/*.coverage.mat']));
  P.summed_cov_track = decell(direc([indir '/*.summed_coverage.fwb']));
end

if iscell(P)
  % multi-set load mode
  if exist('isetname','var'), error('use P.isetname for multi-set load mode'); end
  demand_multiP_files(P);
  M = cell(length(P),1);
  for i=1:length(P)
    M{i} = load_all_mutation_data2(P{i});
  end
  check_multiM_agreement(M);
  return
end

if ~exist('P','var'), P=[]; end

if ~exist('isetname','var'), isetname = 'iset'; end
P=impose_default_value(P,'isetname',isetname);

if isfield(P,'patfile') && ~isfield(P,'patlist')
  P = rename_field(P,'patfile','patlist');
end
if isfield(P,'genefile') && ~isfield(P,'genelist')
  P = rename_field(P,'genefile','genelist');
end
if isfield(P,'mutlist') && ~isfield(P,'mutfile')
  P = rename_field(P,'mutlist','mutfile');
end

P=impose_default_value(P,'mutfile','*required*');
P=impose_default_value(P,'catfile',[]);   % if not specified, will take from covfile
P=impose_default_value(P,'patlist',[]);   % if not specified, will take from covfile; otherwise only field required= "name"
P=impose_default_value(P,'genelist',[]);  % if not specified, will take from covfile; otherwise only field used= "name"
P=impose_default_value(P,'exprcorrfile',[]);          % expression-based BMR corrections, per gene and per category
P=impose_default_value(P,'exprfile',[]);              % expression file, to calculate corrections from
P=impose_default_value(P,'BMR_covariates',[]);        % BMR covariates, for use in neighborhood analysis
P=impose_default_value(P,'exclude_genes_equal_to_zero','none');
P=impose_default_value(P,'ignore_coverage',false);
P=impose_default_value(P,'impute_full_coverage',false);
P=impose_default_value(P,'automatically_add_zero_mutation_patients',true);
if ~P.ignore_coverage
  if ~isfield(P,'covfile') && isfield(P,'terrfile'), P = rename_field(P,'terrfile','covfile'); end
  P=impose_default_value(P,'covfile','*required*');
  P=impose_default_value(P,'ignore_categ_list_mismatch',false);     % NOT RECOMMENDED TO CHANGE THIS!
  tmp = regexprep(P.covfile,'.coverage.mat$','.summed_coverage.fwb');
  P=impose_default_value(P,'summed_cov_track',tmp);      % required for MutSigII features
  P=impose_default_value(P,'patient_min_frac_coverage_required',0.01);
  P=impose_default_value(P,'gene_min_frac_coverage_required',0.25);
else
  P=impose_default_value(P,'summed_cov_track',[]);
end
P=impose_default_value(P,'sample_specific_fwb_directory',[]);
P=impose_default_value(P,'ignore_flank_indels',true);
P=impose_default_value(P,'skip_directly_to_mutsig2_analysis',false);

if ~exist('build','var'), build=[]; end
P=impose_default_value(P,'build',build);

demand_file(P.mutfile);
if ~P.ignore_coverage, demand_file(P.covfile); end
if ~isempty(P.catfile), demand_file(P.catfile); end
if ~isempty(P.patlist), demand_file(P.patlist); end
if ~isempty(P.genelist), demand_file(P.genelist); end
if ~isempty(P.exprcorrfile), demand_file(P.exprcorrfile); end
if ~isempty(P.exprfile), demand_file(P.exprfile); end
if ~isempty(P.summed_cov_track)
  if ~exist(P.summed_cov_track,'file'), fprintf('WARNING:  no summed coverage track available.  MutSig2 will be unhappy!\n'); end
end

if ~isempty(P.exclude_genes_equal_to_zero) && ~strcmpi(P.exclude_genes_equal_to_zero,'none')
  demand_file(P.exclude_genes_equal_to_zero);
end
if ~isempty(P.BMR_covariates) && ~strcmpi(P.BMR_covariates,'none')
  if ischar(P.BMR_covariates), P.BMR_covariates = {P.BMR_covariates}; end
  for i=1:length(P.BMR_covariates), demand_file(P.BMR_covariates{i}); end
end


M = [];
M.name = P.isetname;
M.build = P.build;
M.report = sprintf('LOADING MUTATION DATA ');
if ~isempty(M.name), M.report = [M.report sprintf('FOR SET "%s"',M.name)]; end
M.report = [M.report sprintf('\n')];
fprintf('%s',M.report);

% LOAD ALL DATA

if ~P.ignore_coverage
  % COVERAGE
  demand_file(P.covfile);
  fprintf('Loading coverage... ');
  M.file.cov = P.covfile;
  tmp = load(P.covfile);
  fprintf('\n');
  if ~isfield(tmp,'C') && isfield(tmp,'C1'),  tmp = rename_field(tmp,'C1','C'); end
  if ~isfield(tmp,'C') && isfield(tmp,'D'),  tmp = rename_field(tmp,'D','C'); end
  if ~isfield(tmp,'C'), error('coverage matfile must have either C or C1 or D'); end
  M.cov = tmp.C;

  if ~isfield(M.cov,'sample') || ~(slength(M.cov.sample)>0)
    fprintf('covfile appears to contain no sample-specific information: will use impute_full_coverage\n');
    P.impute_full_coverage = true;
  end

  if P.impute_full_coverage
    demand_field(M.cov,'gene_terr');
  end

  if isfield(M.cov,'gene_totcov') || (isfield(M.cov,'gene_non_cov') && isfield(M.cov,'gene_non_terr'))
    % check for patients with extremely low coverage
    if isfield(M.cov,'gene_non_cov') && isfield(M.cov,'gene_non_terr')   % measure by coding territory if available
      patcov = sum(M.cov.gene_non_cov(:,:,end),1);
      patterr = sum(M.cov.gene_non_terr(:,end),1);
    else  % use total territory
      patcov = sum(M.cov.gene_totcov,1);
      patterr = sum(M.cov.gene.len);
    end
    patfcov = patcov / patterr;
    lowcov = (patfcov<P.patient_min_frac_coverage_required);
    if any(lowcov)
      fprintf('\nThe following %d patients have very low coverage and are being omitted:\n',sum(lowcov));
      disp([as_column(M.cov.sample.name(lowcov)) as_column(num2cell(patfcov(lowcov)))]);
      M.report = [M.report sprintf('Removed %d patients with very low coverage.\n',sum(lowcov))];
      M.cov.sample = reorder_struct(M.cov.sample,~lowcov);
      M.cov.ns = sum(~lowcov);
      if isfield(M.cov,'orig_cov'), M.cov.orig_cov = M.cov.orig_cov(~lowcov,:); end
      if isfield(M.cov,'fcov'), M.cov.fcov = M.cov.fcov(:,~lowcov); end
      if isfield(M.cov,'fcov_coding'), M.cov.fcov_coding = M.cov.fcov_coding(:,~lowcov); end
      if isfield(M.cov,'gene_totcov'), M.cov.gene_totcov = M.cov.gene_totcov(:,~lowcov); end
      if isfield(M.cov,'gene_cov'), M.cov.gene_cov = M.cov.gene_cov(:,~lowcov,:); end
      if isfield(M.cov,'gene_sil_cov'), M.cov.gene_sil_cov = M.cov.gene_sil_cov(:,~lowcov,:); end
      if isfield(M.cov,'gene_non_cov'), M.cov.gene_non_cov = M.cov.gene_non_cov(:,~lowcov,:); end
      if isfield(M.cov,'gene_flank_cov'), M.cov.gene_flank_cov = M.cov.gene_flank_cov(:,~lowcov,:); end
    end
    
    % check for genes with extremely low coverage
    if isfield(M.cov,'gene_non_cov') && isfield(M.cov,'gene_non_terr')   % measure by coding territory if available
      genecov = sum(M.cov.gene_non_cov(:,:,end),2);
      geneterr = M.cov.gene_non_terr(:,end) * M.cov.ns;
    else  % use total territory
      genecov = sum(M.cov.gene_totcov,2);
      geneterr = M.cov.ns * M.cov.gene.len;
    end
    genefcov = genecov ./ geneterr;
    lowcov = (genefcov<P.gene_min_frac_coverage_required);
    if any(lowcov)
      nremove = sum(lowcov);
      fprintf('\nThe following %d genes have very low coverage and are being omitted:\n',nremove);
      disp([as_column(M.cov.gene.name(lowcov)) as_column(num2cell(genefcov(lowcov)))]);
      if nremove>50
        fprintf('The above %d genes have very low coverage and are being omitted.\n',nremove);
      end
      M.report = [M.report sprintf('Removed %d genes with very low coverage.\n',nremove)];
      keeptarg = (~lowcov(M.cov.targ.gidx));
      M.cov.targ = rename_field(M.cov.targ,'gidx','gidx_orig');
      M.cov.targ = reorder_struct(M.cov.targ,keeptarg);
      M.cov.nt = sum(keeptarg);
      if isfield(M.cov,'fcov'), M.cov.fcov = M.cov.fcov(keeptarg,:); end
      if isfield(M.cov,'fcov_coding'), M.cov.fcov_coding = M.cov.fcov_coding(keeptarg,:); end
      M.cov.gene.idx_orig = (1:slength(M.cov.gene))';
      M.cov.gene = reorder_struct(M.cov.gene,~lowcov);
      M.cov.gene.idx_final = (1:slength(M.cov.gene))';
      M.cov.targ.gidx = mapacross(M.cov.targ.gidx_orig,M.cov.gene.idx_orig,M.cov.gene.idx_final);
      M.cov.gene = rmfield(M.cov.gene,{'idx_orig','idx_final'});
      M.cov.ng = sum(~lowcov);
      if isfield(M.cov,'gene_orig_terr'), M.cov.gene_orig_terr = M.cov.gene_orig_terr(~lowcov,:); end
      if isfield(M.cov,'gene_totterr'), M.cov.gene_totterr = M.cov.gene_totterr(~lowcov); end
      if isfield(M.cov,'gene_terr'), M.cov.gene_terr = M.cov.gene_terr(~lowcov,:); end
      if isfield(M.cov,'gene_sil_terr'), M.cov.gene_sil_terr = M.cov.gene_sil_terr(~lowcov,:); end
      if isfield(M.cov,'gene_non_terr'), M.cov.gene_non_terr = M.cov.gene_non_terr(~lowcov,:); end
      if isfield(M.cov,'gene_flank_terr'), M.cov.gene_flank_terr = M.cov.gene_flank_terr(~lowcov,:); end
      if isfield(M.cov,'gene_totcov'), M.cov.gene_totcov = M.cov.gene_totcov(~lowcov,:); end
      if isfield(M.cov,'gene_cov'), M.cov.gene_cov = M.cov.gene_cov(~lowcov,:,:); end
      if isfield(M.cov,'gene_sil_cov'), M.cov.gene_sil_cov = M.cov.gene_sil_cov(~lowcov,:,:); end
      if isfield(M.cov,'gene_non_cov'), M.cov.gene_non_cov = M.cov.gene_non_cov(~lowcov,:,:); end
      if isfield(M.cov,'gene_flank_cov'), M.cov.gene_flank_cov = M.cov.gene_flank_cov(~lowcov,:,:); end
    end
  end
end

if ~isempty(P.summed_cov_track)
  M.file.summed_cov_track = P.summed_cov_track;
end

% GENES
if ~isempty(P.genelist)
  fprintf('Loading gene list... ');
  M.file.gene = P.genelist;
  if get_num_cols(P.genelist)==1
    tmp = load_lines(P.genelist);
    if strcmpi(tmp{1},'name'), tmp(1)=[]; end
  else
    tmp = load_struct(P.genelist);
    tmp = tmp.name;
  end
  fprintf('\n');
  M.gene.name = unique(tmp);
  M.gene.name = M.gene.name(~cellfun('isempty',M.gene.name));

  if ~P.ignore_coverage
    if isfield(M.cov,'gene')
      M.gene.cov_gidx = listmap(M.gene.name,M.cov.gene.name);
      missing = (isnan(M.gene.cov_gidx));
    else
      M.cov.targ.gidx = listmap(M.cov.targ.gene,M.gene.name);
      h = histc(M.cov.targ.gidx,1:slength(M.gene));
      missing = (h==0);
    end
    if any(missing)
      fprintf('\nThe following %d genes are being omitted because coverage is zero/redacted:\n',sum(missing));
      disp(M.gene.name(missing));
      M.gene = reorder_struct(M.gene,~missing);
      M.cov.targ.gidx = listmap(M.cov.targ.gene,M.gene.name);
    end
    % now remove unused targets and genes from M.cov
    if any(isnan(M.gene.cov_gidx)), fprintf('What?\n'); keyboard; end
    gidx = M.gene.cov_gidx;
    % TODO: replace this copy-paste code (and the analogous code for patients) with a subfunction
    if isfield(M.cov,'gene') M.cov.gene = reorder_struct(M.cov.gene,gidx); end
    M.cov.ng = length(gidx);
    if isfield(M.cov,'gene_orig_terr') M.cov.gene_orig_terr = M.cov.gene_orig_terr(gidx,:); end
    if isfield(M.cov,'gene_totterr') M.cov.gene_totterr = M.cov.gene_totterr(gidx,:); end
    if isfield(M.cov,'gene_terr') M.cov.gene_terr = M.cov.gene_terr(gidx,:); end
    if isfield(M.cov,'gene_sil_terr') M.cov.gene_sil_terr = M.cov.gene_sil_terr(gidx,:); end
    if isfield(M.cov,'gene_non_terr') M.cov.gene_non_terr = M.cov.gene_non_terr(gidx,:); end
    if isfield(M.cov,'gene_flank_terr'), M.cov.gene_flank_terr = M.cov.gene_flank_terr(gidx,:); end
    if isfield(M.cov,'gene_totcov'), M.cov.gene_totcov = M.cov.gene_totcov(gidx,:); end
    if isfield(M.cov,'gene_cov'), M.cov.gene_cov = M.cov.gene_cov(gidx,:,:); end
    if isfield(M.cov,'gene_sil_cov'), M.cov.gene_sil_cov = M.cov.gene_sil_cov(gidx,:,:); end
    if isfield(M.cov,'gene_non_cov'), M.cov.gene_non_cov = M.cov.gene_non_cov(gidx,:,:); end
    if isfield(M.cov,'gene_flank_cov'), M.cov.gene_flank_cov = M.cov.gene_flank_cov(gidx,:,:); end
    if isfield(M.cov,'targ')
      M.cov.targ.gidx = listmap(M.cov.targ.gene,M.cov.gene.name);
      M.cov.targ = reorder_struct(M.cov.targ,~isnan(M.cov.targ.gidx));
    end    
  end
else % no genelist provided: use unaltered genelist from coverage file
  if ~P.ignore_coverage
    if ~isfield(M.cov,'gene'), error('cov needs gene field to use this option'); end
    M.gene = [];
    M.gene.name = M.cov.gene.name;
    M.gene.cov_gidx = (1:slength(M.cov.gene))';
    fprintf('Using %d genes in coverage file.\n',slength(M.gene));
  else
    error('Need a genelist if we''re ignoring coverage');
  end
end
M.gene.longname = get_longnames(M.gene.name);

% PATIENTS
if ~isempty(P.patlist)
  fprintf('Loading patient list... ');
  M.file.patient = P.patlist;
  if get_num_cols(P.patlist)==1
    tmp = load_lines(P.patlist);
    if strcmpi(tmp{1},'name'), tmp(1)=[]; end
    tmp2 = []; tmp2.name = tmp;
  else
    tmp2 = load_struct(P.patlist);
  end
  fprintf('\n');
  [u ui] = unique(tmp2.name);
  tmp2 = reorder_struct(tmp2,ui);
  M.patient = tmp2;
  if ~P.ignore_coverage && ~P.impute_full_coverage
    M.cov.sample.pidx = listmap(M.cov.sample.name,M.patient.name);
    h = histc(M.cov.sample.pidx,1:slength(M.patient));
    if any(h>1)
      fprintf('\nThe following patients have multiple entries in coverage table:\n');
      disp(M.patient.name(h>1));
      error('Need to update code to account for this possibility');
    end
    missing = (h==0);
    if any(missing)
      fprintf('\nThe following patients have no coverage information and are being omitted:\n');
      disp(M.patient.name(missing));
      M.report = [M.report sprintf('Removed %d patients with missing coverage information.\n',sum(missing))];
      M.patient = reorder_struct(M.patient,~missing);
      M.cov.sample.pidx = listmap(M.cov.sample.name,M.patient.name);
    end
    M.patient.cidx = listmap(M.patient.name,M.cov.sample.name);
    if P.impute_full_coverage
      M.patient.totcov =  repmat(sum(M.cov.gene_totterr),length(M.patient.cidx),1);
    else
      if isfield(M.cov,'gene_totcov')
        M.patient.totcov = sum(M.cov.gene_totcov(:,M.patient.cidx),1)';
      elseif isfield(M.cov,'totcov')
        M.patient.totcov = sum(M.cov.totcov(:,M.patient.cidx),1)';
      else
        error('Can''t find field for total coverage in M.cov.  Try setting P.impute_full_coverage = true');
      end
    end
    M.patient.cidx = listmap(M.patient.name,M.cov.sample.name);
    M.patient = rmfield(M.patient,'name');
    M.patient = merge_structs({M.patient,reorder_struct(M.cov.sample,M.patient.cidx)});
    % remove coverage information from M.cov for patients outside the patient set
    M.cov.sample = reorder_struct(M.cov.sample,M.patient.cidx);
    M.cov.ns = slength(M.cov.sample);
    if isfield(M.cov,'orig_cov'), M.cov.orig_cov = M.cov.orig_cov(M.patient.cidx,:); end
    if isfield(M.cov,'fcov'), M.cov.fcov = M.cov.fcov(:,M.patient.cidx); end
    if isfield(M.cov,'fcov_coding'), M.cov.fcov_coding = M.cov.fcov_coding(:,M.patient.cidx); end
    if isfield(M.cov,'gene_totcov'), M.cov.gene_totcov = M.cov.gene_totcov(:,M.patient.cidx); end
    if isfield(M.cov,'gene_cov'), M.cov.gene_cov = M.cov.gene_cov(:,M.patient.cidx,:); end
    if isfield(M.cov,'cov'), M.cov.cov = M.cov.cov(:,M.patient.cidx,:); end
    if isfield(M.cov,'totcov'), M.cov.totcov = M.cov.totcov(:,M.patient.cidx); end
    M.patient.cidx = (1:M.cov.ns)';
  end
else % no patient list provided: use unaltered patient list LATER from mutation list and/or coverage file
  fprintf('No patient list provided: will use full list of patients\n');
end

% CATEGORIES
if ~isempty(P.catfile)
  fprintf('Loading categories... ');
  M.file.categ = P.catfile;
  M.categ = load_struct(P.catfile);
  if isfield(M.categ,'rate'), M.categ = make_numeric(M.categ,'rate'); end
  fprintf('\n');
else % no category list provided: use unaltered list from coverage file
  if ~P.ignore_coverage
    M.categ = M.cov.categ;
    fprintf('Using category list from coverage file\n');
  else
    error('Need category file if we''re ignoring coverage.');
  end
end
% legacy cases
idx = grepi('indel|null',M.categ.type,1);
if ~isempty(idx), M.categ.type(idx) = repmat({'non-point'},length(idx),1); end
if isfield(M.categ,'name')
  idx = grepi('indel|null',M.categ.name,1);
  if ~isempty(idx), M.categ.type(idx) = repmat({'non-point'},length(idx),1); end
end
M.categ.typeidx = listmap(lower(M.categ.type),{'point','non-point'});
if any(isnan(M.categ.typeidx))
  error('Unknown category type: %s\n',M.categ.type{find(isnan(M.categ.typeidx),1)});
end
if ~issorted(M.categ.typeidx)
  error('Non-point (indel/null) categories need to be last in the list');
end
if ~P.ignore_coverage
  % check agreement between coverage categories and category list for mutation assignment
  if ~isfield(M.cov,'categ')
    if M.cov.ncat ~= slength(M.categ)
      error('covfile and catfile must have same number of categories');
    end
    fprintf('WARNING: covfile lacks categ field; cannot rigorously verify consistency of categories\n');
  else
    % legacy cases
    idx = grepi('indel|null',M.cov.categ.type,1);
    if ~isempty(idx), M.cov.categ.type(idx) = repmat({'non-point'},length(idx),1); end
    if isfield(M.cov.categ,'name')
      idx = grepi('indel|null',M.cov.categ.name,1);
      if ~isempty(idx), M.cov.categ.type(idx) = repmat({'non-point'},length(idx),1); end
    end
    if ~P.ignore_categ_list_mismatch
      if slength(M.categ)~=slength(M.cov.categ) || ...
            ~all(strcmp(M.categ.name,M.cov.categ.name)) || ...
            ~all(strcmp(M.categ.type,M.cov.categ.type))
        fprintf('\nERROR: covfile and catfile have different categories:\n');
        fprintf('\ncovfile:\n'); look(M.cov.categ);
        fprintf('\ncatfile:\n'); look(M.categ);
        error('Please resolve this discrepancy and try again.');
      end
    end
  end
end


% EXPRESSION CORRECTION FACTORS with a pre-calculated factors file
%    this is the old way, by using an already made .mat exprcorr file 
if ~isempty(P.exprcorrfile)
  if ~isempty(P.exprfile)
    error(['Error: either an input file of expression factors can be a supplied,'...
          'or they can be calculated internally; not both.']);
  end
  fprintf('Loading expression-based BMR correction factors... ');
  load(P.exprcorrfile,'E'); require_fields(E,{'gene','effect'});
  if size(E.effect,2)>1
    % it's multi-column:  check to make sure the categs match
    require_field(E,'categ');
    if size(E.effect,2)~=slength(E.categ)
      error('size(E.effect,2)~=slength(E.categ)');
    end
    if slength(E.categ)~=slength(M.categ) || ~all(strcmpi(M.categ.name,E.categ.name))
      fprintf('ERROR: catfile and exprcorrfile have different categories:\n');
      fprintf('\ncatfile:\n'); look(M.categ);
      fprintf('\nexprcorrfile:\n'); look(E.categ);
      error('Please resolve this discrepancy and try again.');
    end
  end
  default_value = 1;
  M.exprcorr = mapacross(M.gene.name,E.gene.name,E.effect,default_value);
  fprintf('\n');
end

% BMR COVARIATES
% for use in neighborhood analysis
if ~isempty(P.BMR_covariates)
  fprintf('Loading BMR covariate info.\n');
  covar = P.BMR_covariates;
  if ischar(covar), covar = { covar }; end
  M.V = [];
  M.V.file = covar;
  M.V.name = cell(length(covar),1);
  M.V.val = cell(length(covar),1);
  for i=1:length(covar)
    fname = covar{i};
    [z M.V.name{i} z z] = fileparts(fname);
    tmp = load_struct_noheader(fname,{'name','val'});
    if strcmpi(tmp.name{1},'name') || strcmpi(tmp.name{1},'gene'), tmp = reorder_struct_exclude (tmp,1); end
    tmp = make_numeric(tmp,'val');
    M.V.val{i} = nansub(tmp.val,listmap(M.gene.name,tmp.name));
  end
end

% GENES TO IGNORE IF ZERO (RPKM or other metric)
if ~isempty(P.exclude_genes_equal_to_zero) & ~strcmp(P.exclude_genes_equal_to_zero, 'none')
  fprintf('Loading exclude_genes_equal_to_zero file.\n');
  fname = P.exclude_genes_equal_to_zero;
  M.file.exclude_genes_equal_to_zero = fname;
  tmp = load_struct_noheader(fname,{'name','val'});
  if strcmpi(tmp.name{1},'name') || strcmpi(tmp.name{1},'gene'), tmp = reorder_struct_exclude (tmp,1); end
  tmp = make_numeric(tmp,'val');
  tmp = reorder_struct(tmp,tmp.val==0);
  gidx = listmap(tmp.name,M.gene.name); gidx(isnan(gidx)) = [];
  M.gene.exclude_because_zero = false(slength(M.gene),1);
  M.gene.exclude_because_zero(gidx) = true;
end

% MUTATIONS
fprintf('Loading mutations... ');
M.file.mut = P.mutfile;
M.mut = load_struct(P.mutfile);
fprintf('\n');
M.report = sprintf('Number of mutations before filtering:\t\t%d\n',slength(M.mut));

fprintf('Processing mutation list...\n');
M.mut = add_and_convert_simple_fieldnames(M.mut,P);
M.mut = add_helper_is_fields(M.mut);

if isfield(M.mut,'firehose_patient_id'), M.mut.patient = M.mut.firehose_patient_id; end   % MAKE SURE

% make sure we know patient name for each mutation
if ~isfield(M.mut,'patient')
  error('Can''t find patient names in the MAF!  Please provide in either "patient" or "Tumor_Sample_Barcode"');
end

% if no patient list was provided, extract one now from the mutation list
if ~isfield(M,'patient')
  pm = unique(M.mut.patient);
  if isfield(M.cov,'sample') && isfield(M.cov.sample,'name') && slength(M.cov.sample)>0
    pc = M.cov.sample.name;
    if length(intersect(pc,pm))==0
      fprintf('Patient names in covfile don''t match patient names in maffile: will use maffile.\n');
    else
      if P.automatically_add_zero_mutation_patients
        blankpat = setdiff(pc,pm);
        if ~isempty(blankpat)
          fprintf(['Covfile contains %d patients not in mutation list:'...
                   'adding these "zero-mutations" patients to analysis:\n'],length(blankpat));
          disp(blankpat);
          pm = union(pm,blankpat);
        end
      end
    end
  end
  M.patient = [];
  M.patient.name = pm;
  M.patient.cidx = listmap(M.patient.name,M.cov.sample.name);
end

% make sure the mutations are adequately mapped to the patient list
z = listmap(M.mut.patient,M.patient.name);
if all(isnan(z))
  % stopgap measure: try truncating the TCGA names
  M.mut.patient = regexprep(M.mut.patient,'^(TCGA-..-....).*$','$1');
  M.patient.name = regexprep(M.patient.name,'^(TCGA-..-....).*$','$1');
  z = listmap(M.mut.patient,M.patient.name);
end

if all(isnan(z))
  error('Couldn''t find patient names in the MAF!');
end
if mean(isnan(z))>0.1
  fprintf('*** WARNING ****   only %d%% of mutations were mapped to the patient list\n',100*mean(~isnan(z)));
end

% add other necessary fields to M.mut (if they weren't already added by MutSigPreprocess)
if ~isfield(M.mut,'context65')
  fprintf('context65 was missing from MAF:  loading it now...\n');
  if isfield(P,'contextdir')
    fprintf('Using P.contextdir\n');
    contextdir = P.contextdir;
  elseif contains(P.build,'hg18')
    fprintf('Using hardcoded default path: will not work on external systems.\n');
    contextdir = '/cga/tcga-gsc/home/lawrence/db/hg18/context65';
  elseif contains(P.build,'hg19')
    fprintf('Using hardcoded default path: will not work on external systems.\n');
    contextdir = '/cga/tcga-gsc/home/lawrence/db/hg19/context65';
  else
    error('no context65, and unknown build');
  end
  M.mut.context65 = get_context(M.mut.chr,M.mut.start,contextdir);
end

if wouldbenumeric(M.mut.gene) && isfield(M.mut,'gene_name')
  M.mut.gene = M.mut.gene_name;
end
if wouldbenumeric(M.mut.patient) && isfield(M.mut,'patient_name')
  M.mut.patient = M.mut.patient_name;
end

% remove any "non-mutations"
idx = find(strcmpi(M.mut.newbase,M.mut.ref_allele));
if ~isempty(idx)
  fprintf('Removing the following %d non-mutations:\n',length(idx));
  look(M.mut,idx);
  M.mut = reorder_struct_exclude(M.mut,idx);
  M.report = [M.report sprintf('After removing %d non-mutations:\t%d\n',length(idx),slength(M.mut))];
end

% assign mutations to categories
[M.mut.categ M.mut.categ_ignoring_null_categ] = assign_mut_categs(M.mut,M.categ);
double_null_categ = grepi('double.?null',M.categ.name,1);
if ~isempty(double_null_categ)
  num_double_null = sum(M.mut.categ==double_null_categ(1));
else
  num_double_null = 0;
end

% choose mutations that are in the patient/gene/category set

M.ng = slength(M.gene);
M.np = slength(M.patient);
M.mut.gene = regexprep(M.mut.gene,'^([^\|]+).*','$1');
M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);
M.mut.pat_idx = listmap(M.mut.patient,M.patient.name);
M.mut.categ(isnan(M.mut.categ) | M.mut.categ<1 | M.mut.categ>slength(M.categ)) = nan;
M.mut.categ(isnan(M.mut.categ_ignoring_null_categ) | M.mut.categ_ignoring_null_categ<1 | ...
                                 M.mut.categ_ignoring_null_categ>slength(M.categ)) = nan;
M.use = find(~isnan(M.mut.gene_idx) & ~isnan(M.mut.pat_idx) & ~isnan(M.mut.categ));
pidx = find(isnan(M.mut.pat_idx));
gidx = find(isnan(M.mut.gene_idx));
cidx = find(isnan(M.mut.categ));
if ~isempty(pidx)
  fprintf('Patients not in patient set:');
  count(M.mut.patient(pidx),1);
  fprintf('Mutations not in patient set:');
  if length(unique(M.mut.gene(pidx)))>100, fprintf('   [too many to list]\n'); else count(M.mut.gene(pidx),1); end
end
if ~isempty(gidx)
  fprintf('Mutations not in gene set:');
  if length(unique(M.mut.gene(gidx)))>100, fprintf('   [too many to list]\n'); else count(M.mut.gene(gidx),1); end
end
if ~isempty(cidx)
  fprintf('Categories not in categ set:');
  count(M.mut.categ(cidx),1);
  fprintf('Mutations not in categ set:');
  if length(unique(M.mut.gene(cidx)))>100, fprintf('   [too many to list]\n'); else count(M.mut.gene(cidx),1); end
end
fprintf('\nMutations in this patient+gene+categ set: %d/%d\n', length(M.use),slength(M.mut));
fprintf('\t%d mutations not in patient set\n',length(pidx));
fprintf('\t%d mutations not in gene set\n',length(gidx));
fprintf('\t%d mutations not in categ set\n',length(cidx));
if num_double_null>0
  fprintf('\t\t(Note: at least %d of these are excluded because they are partners to double-nulls.)\n',num_double_null);
end
f = length(cidx) / slength(M.mut);
if f>0.05, fprintf('WARNING: %d%% of mutations are outside the category set!\n',round(100*f)); end
if isempty(M.use), fprintf('No mutations!'); keyboard; error('No mutations!'); end

nxp = length(pidx);
nxg = length(setdiff(gidx,pidx));
nxc = length(setdiff(cidx,[pidx;gidx]));
if nxp>0, M.report = [M.report sprintf('After removing %d mutations outside patient set:\t%d\n',nxp,slength(M.mut)-nxp)]; end
if nxg>0, M.report = [M.report sprintf('After removing %d mutations outside gene set:\t%d\n',nxg,slength(M.mut)-nxp-nxg)]; end
if nxc>0, M.report = [M.report sprintf('After removing %d mutations outside category set:\t%d\n',nxc,slength(M.mut)-nxp-nxg-nxc)]; end 

if any(strcmpi(M.categ.name,'Total')), error('catfile should not include entry for "Total"'); end
M.mutclass = [M.categ.name' 'Total'];
M.TOT = length(M.mutclass);
M.NUM_INDEL_CLASSES = sum(M.categ.typeidx==2);

if P.ignore_flank_indels
  idx = M.use(M.mut.is_flank(M.use) & M.mut.is_indel(M.use));
  fprintf('Ignoring %d flank indels\n', length(idx));
  M.use = setdiff(M.use,idx);
end

% get rid of the mutations we're not going to use
M.mut = reorder_struct(M.mut,M.use);
M.use = as_row(1:slength(M.mut));

% if P.sample_specific_fwb_directory was provided, add this to M.patient
if ~isempty(P.sample_specific_fwb_directory)
  fprintf('Pointing to sample-specific FWBs in directory %s\n',P.sample_specific_fwb_directory);
  M.patient.fwb = regexprep(M.patient.name,'^(.*)$',[P.sample_specific_fwb_directory '/$1.cov.fwb']);
end

% if P.mutsig2_power_calculation_in_dir was provided, add this to M.patient 
if isfield(P,'mutsig2_power_calculation_in_dir') && ~isempty(P.mutsig2_power_calculation_in_dir)
  fprintf('Pointing to sample-specific exact read-count directory %s\n', P.mutsig2_power_calculation_in_dir); 
  M.patient.power_files =  regexprep(M.patient.name,'^(.*)$',[P.mutsig2_power_calculation_in_dir '/$1-Tumor.cov.fwb']);
  fprintf('Verifying that files exist... \n');
  ok = demand_files(M.patient.power_files); 
  fprintf('%d out of %d .fwb files were found...', sum(ok), length(M.patient.power_files));
  if sum(ok) == 0
    fprintf('No .fwb files were found; will check for .cbb files...\n');
    M.patient.power_files = regexprep(M.patient.name,'^(.*)$',[P.mutsig2_power_calculation_in_dir '/$1-Tumor.cov.cbb']);
    ok = demand_files(M.patient.power_files);
    fprintf('%d out of %d .cbb files were found...', sum(ok), length(M.patient.power_files));
  end 
  if sum(ok) == 0 
    error('No files were found in the provided read count directory, power calculation mode cannot be performed');
  end 

end 

if P.skip_directly_to_mutsig2_analysis
  return
end


% build territory tables (N_terr)
if isfield(M.cov,'gene_terr'), M.N_terr = M.cov.gene_terr(M.gene.cov_gidx,:); end
if isfield(M.cov,'gene_sil_terr'), M.N_sil_terr = M.cov.gene_sil_terr(M.gene.cov_gidx,:); end
if isfield(M.cov,'gene_non_terr'), M.N_non_terr = M.cov.gene_non_terr(M.gene.cov_gidx,:); end
if isfield(M.cov,'gene_flank_terr'), M.N_flank_terr = M.cov.gene_flank_terr(M.gene.cov_gidx,:); end
if M.NUM_INDEL_CLASSES>0   % (can use totals from indel class)
  idx = M.TOT-M.NUM_INDEL_CLASSES;
  for c=M.TOT-M.NUM_INDEL_CLASSES+1:M.TOT
    if isfield(M,'N_terr'), M.N_terr(:,c) = M.N_terr(:,idx); end
    if isfield(M,'N_sil_terr'), M.N_sil_terr(:,c) = M.N_sil_terr(:,idx); end
    if isfield(M,'N_non_terr'), M.N_non_terr(:,c) = M.N_non_terr(:,idx); end
    if isfield(M,'N_flank_terr'), M.N_flank_terr(:,c) = M.N_flank_terr(:,idx); end
  end
end

% build coverage tables (N_cov)
if ~(P.ignore_coverage || P.impute_full_coverage)

  fprintf('\nBuilding coverage table...');
  
  M.N_cov = zeros(M.ng,M.TOT,M.np);

  if isfield(M.cov,'gene_sil_cov'), M.N_sil_cov = zeros(M.ng,M.cov.ncat,M.np); end
  if isfield(M.cov,'gene_non_cov'), M.N_non_cov = zeros(M.ng,M.cov.ncat,M.np); end
  if isfield(M.cov,'gene_flank_cov'), M.N_flank_cov = zeros(M.ng,M.cov.ncat,M.np); end

  if ~isfield(M.cov,'gene')
    % old method (coverage table doesn't include gene-collapsed version)
    for g=1:M.ng
      idx = find(M.cov.targ.gidx==g);
      M.gene.chr(g,1) = M.cov.targ.chr(idx(1));
      M.gene.start(g,1) = min(M.cov.targ.start(idx));
      M.gene.end(g,1) = max(M.cov.targ.end(idx));
      M.gene.len(g,1) = sum(M.cov.targ.len(idx));
      M.gene.gc(g,1) = weighted_mean(M.cov.targ.gc(idx),M.cov.targ.len(idx));
      if isfield(M.cov.targ,'type'), M.gene.type(g,1) = M.cov.targ.type(idx(1)); end
      for c=1:M.cov.ncat, M.N_cov(g,c,:) = sum(M.cov.cov(idx,M.patient.cidx,c),1); end
      for c=M.TOT-M.NUM_INDEL_CLASSES:M.TOT, M.N_cov(g,c,:) = sum(M.cov.totcov(idx,M.patient.cidx),1); end
    end

  else  % new method (coverage table already includes gene-collapsed version)

    M.gene = merge_structs({M.gene,reorder_struct(M.cov.gene,M.gene.cov_gidx)});
    for c=1:M.cov.ncat, M.N_cov(:,c,:) = M.cov.gene_cov(M.gene.cov_gidx,M.patient.cidx,c); end
    for c=M.TOT-M.NUM_INDEL_CLASSES:M.TOT, M.N_cov(:,c,:) = M.cov.gene_totcov(M.gene.cov_gidx,M.patient.cidx); end

    for c=1:M.cov.ncat     
      if isfield(M.cov,'gene_sil_cov'), M.N_sil_cov(:,c,:) = M.cov.gene_sil_cov(M.gene.cov_gidx,M.patient.cidx,c); end
      if isfield(M.cov,'gene_non_cov'), M.N_non_cov(:,c,:) = M.cov.gene_non_cov(M.gene.cov_gidx,M.patient.cidx,c); end
      if isfield(M.cov,'gene_flank_cov'), M.N_flank_cov(:,c,:) = M.cov.gene_flank_cov(M.gene.cov_gidx,M.patient.cidx,c); end
    end
    if M.NUM_INDEL_CLASSES>0   % (can use totals from indel class)
      idx = M.TOT-M.NUM_INDEL_CLASSES;
      for c=M.TOT-M.NUM_INDEL_CLASSES+1:M.TOT
        if isfield(M.cov,'gene_sil_cov'), M.N_sil_cov(:,c,:) = M.N_sil_cov(:,idx,:); end
        if isfield(M.cov,'gene_non_cov'), M.N_non_cov(:,c,:) = M.N_non_cov(:,idx,:); end
        if isfield(M.cov,'gene_flank_cov'), M.N_flank_cov(:,c,:) = M.N_flank_cov(:,idx,:); end
      end
    else
      error('Can''t compute totals: no indel-class coverage available!');
    end
  end
end

% build mutation tables (n)

fprintf('\nBuilding mutation tables\n');

M.use_silent = [];
M.use_nonsilent = [];
imposs = [];

if any(M.mut.is_flank)
  M.use_flank = [];
  M.n_flank = zeros(M.ng,M.TOT,M.np);
  use_flank = true;
else
  use_flank = false;
end

M.n_silent = zeros(M.ng,M.TOT,M.np);
M.n_nonsilent = zeros(M.ng,M.TOT,M.np);
M.n_nonsilent_ignoring_null_categ = zeros(M.ng,M.TOT,M.np);

method=2;
if method==2
  % new, fast method
  M.use_nonsilent = as_row(M.use(M.mut.is_coding(M.use) & ~M.mut.is_silent(M.use)));
  M.use_silent = as_row(M.use(M.mut.is_coding(M.use) & M.mut.is_silent(M.use)));
  if use_flank
    M.use_flank = as_row(M.use(M.mut.is_flank(M.use)));
  end
  for c=1:slength(M.categ)
    idx = M.use_nonsilent(M.mut.categ(M.use_nonsilent)==c);
    M.n_nonsilent(:,c,:) = hist2d_fast(M.mut.gene_idx(idx),M.mut.pat_idx(idx),1,M.ng,1,M.np);
    idx = M.use_silent(M.mut.categ(M.use_silent)==c);
    M.n_silent(:,c,:) = hist2d_fast(M.mut.gene_idx(idx),M.mut.pat_idx(idx),1,M.ng,1,M.np);
    idx = M.use_nonsilent(M.mut.categ_ignoring_null_categ(M.use_nonsilent)==c);
    M.n_nonsilent_ignoring_null_categ(:,c,:) = hist2d_fast(M.mut.gene_idx(idx),M.mut.pat_idx(idx),1,M.ng,1,M.np);
    if use_flank
      idx = M.use_flank(M.mut.categ(M.use_flank)==c);
      M.n_flank(:,c,:) = hist2d_fast(M.mut.gene_idx(idx),M.mut.pat_idx(idx),1,M.ng,1,M.np);
    end
  end
  s = M.n_nonsilent+M.n_silent+M.n_nonsilent_ignoring_null_categ;
  if isfield(M,'n_flank'), s = s + M.n_flank; end
  if isfield(M,'N_cov')
    imposs = find(M.N_cov==0 & s>0);
    if ~isempty(imposs)
      M.n_nonsilent(imposs)=0; M.n_nonsilent_ignoring_null_categ(imposs)=0;
      M.n_silent(imposs)=0;
      if isfield(M,'n_flank'), M.n_flank(imposs)=0; end
      fprintf('Omitting %d "impossible" mutations in gene-patient-category bins having zero coverage.\n',length(imposs));
      M.report = [M.report ...
                  sprintf('After removing %d "impossible" mutations in\ngene-patient-category bins of zero coverage:\t%d\n',...
                          length(imposs),slength(M.mut)-nxc-nxp-nxg-length(imposs))];
      imposs=[];
    end
  end

elseif method==1
  % old, slow method
  for i=1:length(M.use), if ~mod(i,10000), fprintf('%d/%d ',i,length(M.use)); end
    m = M.use(i);
    if P.ignore_flank_indels && M.mut.is_flank(m) && M.mut.is_indel(m), continue; end
    p = M.mut.pat_idx(m);
    g = M.mut.gene_idx(m);
    c = M.mut.categ(m);
    c_ignoring_null_categ = M.mut.categ_ignoring_null_categ(m);
    if isnan(p) || isnan(g), error('"use" includes mutations outside patient/gene set'); end
    if c<1 || c>slength(M.categ), error('"use" includes mutations outside category set'); end
    if ~(P.ignore_coverage || P.impute_full_coverage)
      if M.N_cov(g,c,p)==0 || M.N_cov(g,c_ignoring_null_categ,p)==0   % 'impossible mutation'
        imposs(end+1) = m;
        continue;
    end,end
    if M.mut.is_flank(m)
      M.use_flank(end+1) = m;
      M.n_flank(g,c,p) = M.n_flank(g,c,p) + 1;
    elseif M.mut.is_coding(m)
      if M.mut.is_silent(m)
        M.use_silent(end+1) = m;
        M.n_silent(g,c,p) = M.n_silent(g,c,p) + 1;
      else % nonsilent
        M.use_nonsilent(end+1) = m;
        M.n_nonsilent(g,c,p) = M.n_nonsilent(g,c,p) + 1;
        M.n_nonsilent_ignoring_null_categ(g,c_ignoring_null_categ,p) = ...
            M.n_nonsilent_ignoring_null_categ(g,c_ignoring_null_categ,p) + 1;
  end,end,end, if length(M.use)>=10000, fprintf('\n'); end

end % method 

if isfield(M,'n_flank'), M.n_flank(:,M.TOT,:) = sum(M.n_flank,2); end
M.n_silent(:,M.TOT,:) = sum(M.n_silent,2);
M.n_nonsilent(:,M.TOT,:) = sum(M.n_nonsilent,2);
M.n_nonsilent_ignoring_null_categ(:,M.TOT,:) = sum(M.n_nonsilent_ignoring_null_categ,2);

if isfield(M,'use_flank')
  M.use = unique([M.use_flank,M.use_silent,M.use_nonsilent]);
else
  M.use = unique([M.use_silent,M.use_nonsilent]);
end

if ~isempty(imposs)    % remove impossible mutations
  fprintf('Omitting %d "impossible" mutations in gene-patient-category bins having zero coverage:\n',length(imposs));
  count(M.mut.gene(imposs),1);
  M.report = [M.report ...
    sprintf('After removing %d "impossible" mutations in\ngene-patient-category bins of zero coverage:\t%d\n',...
      length(imposs),slength(M.mut)-nxc-nxp-nxg-length(imposs))];
end


% CALCULATE EXPRESSION CORRECTION FACTORS // added on 11/20/2011 by Petar, to calculate the expression correlation factors in place  
if ~isempty(P.exprfile)
  if ~isempty(P.exprcorrfile)
    error('Error: either an input file of expression factors can be a supplied, or they can be calculated internally; not both.');
  end
  fprintf('Loading expression info.\n');

  %% Preprocess mutation and expression data for factor calculation 
  [cov_matrix_tot_nopats_nocats, mut_matrix_tot_nopats_nocats, E1, M2] =  preprocess_expression(M, P.exprfile);
  G = struct();
  G.nfit = mut_matrix_tot_nopats_nocats;
  G.Nfit = cov_matrix_tot_nopats_nocats;
  V = struct();
  V.val = E1.expr;
  
  G = analyze_mutrate_covariates(G, V);

  %% Build exprcorr factors
  
  EC = struct();
  EC.categ = M.categ;
  EC.gene = M2;
  EC.effect = repmat(G.Ffit, 1, slength(M.categ));
  

  if slength(EC.categ)~=slength(M.categ) || ~all(strcmpi(M.categ.name,EC.categ.name))
    fprintf('ERROR: catfile and exprcorrfile have different categories:\n');
    fprintf('\ncatfile:\n'); look(M.categ);
    fprintf('\nexprcorrfile:\n'); look(EC.categ);
    error('Please resolve this discrepancy and try again.');
  end

  M.exprcorr = ones(slength(M.gene),slength(M.categ));
  idx = listmap(M.gene.name,EC.gene.name);
  idx2 = find(~isnan(idx));
  M.exprcorr(idx2,:) = EC.effect(idx(idx2),:);
  fprintf('\n');
end



fprintf('\nAll data loaded\n');
