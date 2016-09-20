function perform_COSMIC_overlap_analysis(M,P,COS)
% perform_COSMIC_overlap_analysis(M[,P,COS])
%
% Mike Lawrence 2008-2011

if ~exist('P','var'),P=[]; end
P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'build_cosmic',P.build);
P = impose_default_value(P,'patients_to_include','all');
P = impose_default_value(P,'cosmic_file',[]);
P = impose_default_value(P,'match_margin',2);
P = impose_default_value(P,'include_silent_mutations',false);
P = impose_default_value(P,'min_reports_in_cosmic_to_consider_mutation',1);
P = impose_default_value(P,'cosmic_report_filename','*required*');
P = impose_default_value(P,'cosmic_report2_filename','*required*');
P = impose_default_value(P,'cosmic_mutations_outname','*required*');

if ~exist('COS','var'), COS = load_cosmic_database(P); end

% n = number of mutations at COSMIC sites
% N = number of unique COSMIC sites x number of patients
% rate = from mut_analysis
demand_field(M,'mutrate');
if isfield(M.mutrate,'tot')
  rate = M.mutrate.tot.hat;
elseif isfield(M.mutrate,'ss')
  fprintf('ssBMR''s not yet implemented in COSMIC analysis: using the mean rate\n');
  rate = mean(M.mutrate.ss.tot.hat);
else
  error('Can''t find mutation rates for COSMIC analysis!');
end

% choose eligible mutations

if P.include_silent_mutations
  M.mut = reorder_struct(M.mut,M.use);
else
  COS.mut = reorder_struct(COS.mut,grepvi('silent|synonymous',COS.mut.MutationDescription,1));
  if isfield(M,'use_nonsilent')
    M.mut = reorder_struct(M.mut,M.use_nonsilent);
  elseif isfield(M.mut,'is_coding') && isfield(M.mut,'is_silent')
    M.mut = reorder_struct(M.mut,M.mut.is_coding & ~M.mut.is_silent);
  elseif isfield(M.mut,'type')
    M.mut = reorder_struct_exclude(M.mut,ismember(M.mut.type,{'Synonymous','Silent','Intron','UTR','3''-UTR','5''-UTR'}));
  else
    fprintf('Unable to find information to distinguish nonsilent coding mutations... using all mutations.\n');
  end
end

if ischar(P.patients_to_include) && strcmpi(P.patients_to_include,'all')
  % keep all mutations
else
  if ischar(P.patients_to_include), P.patients_to_include = { P.patients_to_include }; end
  if isnumeric(P.patients_to_include)
    if ~isfield(M.mut.pat_idx)
      M.mut.pat_idx = listmap(M.mut.patient,M.patient.name);
    end
    M.mut = reorder_struct(M.mut,ismember(M.mut.pat_idx,P.patients_to_include));
  elseif iscellstr(P.patients_to_include)
    M.mut = reorder_struct(M.mut,ismember(M.mut.patient,P.patients_to_include));
  else
    error('unknown format for P.patients_to_include');
  end
end

% find COSMIC overlap

M.mut.n_cos = zeros(slength(M.mut),1);

X = [];
X.gene = intersect(COS.gene.name,M.gene.name);
if isfield(M.gene,'longname'), X.description = mapacross(X.gene,M.gene.name,M.gene.longname); end

z  = nan(slength(X),1);
X.n = z;             % number of mutations in this gene in the individual set
X.cos = z;           % number of unique mutated sites in this gene in COSMIC
X.n_cos = z;         % overlap between n and cos
X.N_cos = z;         % number of individuals x cos
X.cos_ev = z;        % cosmic evidence = number of reports in COSMIC at the sites mutated in this gene

cidx = listmap(COS.gene.name(COS.mut.gene),X.gene);
if isnumeric(M.mut.gene)
  midx = listmap(M.gene.name(M.mut.gene),X.gene);
else
  midx = listmap(M.mut.gene,X.gene);
end

fprintf('COSMIC analysis:  gene ');
for g=1:slength(X), if ~mod(g,1000), fprintf('%d/%d ',g,slength(X)); end
  ci = find(cidx==g);
  mi = find(midx==g);
  X.n(g) = length(mi);
  % make list of all basepairs reported at least once in cosmic
  cospos = COS.mut.start(ci);
  idx = find(COS.mut.end(ci)>COS.mut.start(ci));
  for i=1:length(idx)
    cospos = [cospos;(COS.mut.start(ci(idx(i)))+1:COS.mut.end(ci(idx(i))))'];
  end
  cospos = unique(cospos);
  % count how many times each basepair is reported at least once in cosmic
  cospos_ct = zeros(length(cospos),1);
  for i=1:length(ci)
    for j=COS.mut.start(ci(i)):COS.mut.end(ci(i))
      idx = find(cospos==j);
      cospos_ct(idx) = cospos_ct(idx) + 1;
    end
  end
  % impose threshold on how many times a basepair has to be reported in order to count toward this analysis
  threshold = P.min_reports_in_cosmic_to_consider_mutation;
  idx = find(cospos_ct>=threshold);
  cospos = cospos(idx);
  cospos_ct = cospos_ct(idx);
  % finish analysis for this geene
  X.cos(g) = length(cospos);
  X.n_cos(g) = 0;
  for i=1:length(mi)
    idx = intersect(M.mut.start(mi(i))-P.match_margin:M.mut.end(mi(i))+P.match_margin,cospos);
    if ~isempty(idx)
      X.n_cos(g) = X.n_cos(g)+1;
      % also, for this mutation, count how many COSMIC mutations overlap it
      idx = find(cospos-P.match_margin <= M.mut.end(mi(i)) & cospos+P.match_margin >= M.mut.end(mi(i)));
      M.mut.n_cos(mi(i)) = sum(cospos_ct(idx));
    end
  end 
  X.cos_ev(g) = sum(M.mut.n_cos(mi));
end,fprintf('\n');

% first ordering: by p-value for observed amount of overlap
X.N_cos = M.np * X.cos;
X.p = 1-binocdf(X.n_cos-1,X.N_cos,rate);
X.q = calc_fdr_value(X.p);
X.ratio = X.n_cos ./ X.N_cos;
X = sort_struct(X,{'q','p','ratio'},[1 1 -1]);
X = rmfield(X,'ratio');
X.rank = (1:slength(X))';
X = orderfields_first(X,{'rank'});
save_struct(X,P.cosmic_report_filename);

% second ordering: by number of reports in COSMIC
Y = X;
Y = rmfields(Y,{'N_cos','q','p'});
Y = sort_struct(Y,{'cos_ev'},-1);
Y.rank = (1:slength(Y))';
Y = orderfields_first(Y,{'rank'});
save_struct(Y,P.cosmic_report2_filename);

% save list of mutations seen also in COSMIC
O = reorder_struct(M.mut,M.mut.n_cos>0);
O = sort_struct(O,'n_cos',-1);
save_struct(O,[P.cosmic_mutations_outname '_verbose.txt']);

% brief report for Firehose report
if isfield(O,'patient_name'), O.patient = O.patient_name; end
if isfield(O,'gene_name'), O.gene = O.gene_name; end
O = keep_fields_that_exist(O,{'patient','chr','start','end','type',...
                    'gene','Protein_Change','proteinchange','amino_acid_change',...
                    'domain','ucsc_cons','n_cos'});
save_struct(O,P.cosmic_mutations_outname);




