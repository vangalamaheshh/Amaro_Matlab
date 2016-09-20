function M = build_n_and_N_tables(M,P)
% EXTREMELY OBSOLETE
% probably you're looking for expand_to_n_and_N()

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'manually_remove_offending_KRAS_exon',false);
P=impose_default_value(P,'use_partitioned_rollup_coverage',false);

% build mutation tables (n)

fprintf('Building mutation tables\n');

M.mutclass = { ...
    'A' ;...
    'T' ;...
    'C (CpG)' ;...
    'C (TpC)' ;...
    'C (other)' ;...
    'G (CpG)' ;...
    'G (GpA)' ;...
    'G (other)' ;...
    'Inframe_Ins' ;...
    'Frameshift_Ins', ;...
    'Inframe_Del' ;...
    'Frameshift_Del' ;...
    'Frameshift' ;...
    'Total' ;...
};

M.TOT = length(M.mutclass);
M.NUM_INDEL_CLASSES = 5;

M.n_silent = zeros(M.ng,M.TOT,M.np);
M.n_missense = zeros(M.ng,M.TOT,M.np);
M.n_nonsense = zeros(M.ng,M.TOT,M.np);
M.n_splice = zeros(M.ng,M.TOT,M.np);
M.n_indel = zeros(M.ng,M.TOT,M.np);

for i=1:length(M.use)
  m = M.use(i);
  p = M.mut.patient(m);
  g = M.mut.gene(m);
  if isnan(p) || isnan(g), error('"use" includes mutations outside patient/gene set'); end
  if strcmpi(M.mut.class{m},'Indel')
    % INDEL
    class = find(ismember(M.mutclass, M.mut.type{m}));
    if isempty(class)
      fprintf('Mut %d: Unknown indel type %s\n', m, M.mut.type{m});
    else
      M.n_indel(g,class,p) = M.n_indel(g,class,p) + 1;
    end
  else
    % NON-INDEL
    class = find(ismember(M.mutclass, M.mut.class{m}));
    if isempty(class)
      fprintf('Mut %d: Unknown mutation class %s\n', m, M.mut.class{m});
    else
      if strcmpi(M.mut.type{m},'Nonsense')
        M.n_nonsense(g,class,p) = M.n_nonsense(g,class,p) + 1;
      elseif strcmpi(M.mut.type{m},'Synonymous')
        M.n_silent(g,class,p) = M.n_silent(g,class,p) + 1;
      elseif strcmpi(M.mut.type{m},'Missense')
        M.n_missense(g,class,p) = M.n_missense(g,class,p) + 1;
      elseif strncmpi(M.mut.type{m},'Splice',6)
        M.n_splice(g,class,p) = M.n_splice(g,class,p) + 1;
      else
        fprintf('Mut %d: Unknown mutation type %s\n', m, M.mut.type{m});
      end
    end
  end
end

M.n_missense(:,M.TOT,:) = sum(M.n_missense,2);
M.n_silent(:,M.TOT,:) = sum(M.n_silent,2);
M.n_nonsense(:,M.TOT,:) = sum(M.n_nonsense,2);
M.n_splice(:,M.TOT,:) = sum(M.n_splice,2);
M.n_indel(:,M.TOT,:) = sum(M.n_indel,2);
M.n_nonsilent = M.n_missense + M.n_nonsense + M.n_splice + M.n_indel;

% build coverage table (N)

fprintf('Creating coverage tables\n');

% total possible coverage ("territory")

M.N_terr = zeros(M.ng,M.TOT);
tmp = cat(2,M.genecomp.dat{2:end});
tmp = [tmp nan(length(M.genecomp.dat{2}),M.NUM_INDEL_CLASSES) sum(tmp,2)];   % #indels=NaN
for i=1:length(M.genecomp.gene)
  g = M.genecomp.gene(i);
  if ~isnan(g)
    M.N_terr(g,:) = tmp(i,:);
  end
end

% actual coverage

M.N_cov = zeros(M.ng,M.TOT,M.np);

tmp = cat(2,M.cov.dat{3:end});
tmp = [tmp nan(length(M.cov.dat{3}),M.NUM_INDEL_CLASSES) sum(tmp,2)];    % #indels=NaN
for i=1:length(M.cov.gene)
    g = M.cov.gene(i);
    p = M.cov.patient(i);
    if ~isnan(g) && ~isnan(p)
        M.N_cov(g,:,p) = tmp(i,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          %
% TSP COVERAGE TINKERING   %
%                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MANUALLY REMOVE FIRST EXON OF KRAS FROM GENBOREE DEFINITION

if P.manually_remove_offending_KRAS_exon
  g=grep('KRAS',M.gene.name,1);
  if isempty(g)
    error('Could not locate KRAS in gene list in order to manually remove offending exon!\n');
  else
    fprintf('Manually removing offending first exon of KRAS\n');
    M.N_terr(g,:) = round(M.N_terr(g,:)*(707/900));
    M.N_cov(g,:,:) = round(M.N_cov(g,:,:)*(707/900));
  end
end


% USE GENBOREE-PARTITIONED ROLLUP FILES INSTEAD OF GENBOREE-PROCESSED VERBOSE COVERAGE

if P.use_partitioned_rollup_coverage
  fprintf('Partitioning rollup coverage:\n');
  fprintf('  Total before: %d\n', sum(sum(M.N_cov(:,M.TOT,:))));
  for p=1:M.np
    pidx = find(M.rollup.patient==p);
    if isempty(pidx)
      fprintf('Missing patient %s', M.patient.name(p));
      continue               % should never happen
    end
    for g=1:M.ng
      gidx = find(M.rollup.gene==g);
      if isempty(gidx)
        fprintf('Missing gene %s\n', M.gene.name(g));
        continue
      end       % this problem only happens for NAT6, missing in their list.
      cv = M.rollup.cov(gidx,pidx);
      if M.N_terr(g,M.TOT)>0
        fracs = M.N_terr(g,:)./M.N_terr(g,M.TOT);
      else
        fracs = sum(M.N_terr)./sum(M.N_terr(:,M.TOT));
      end
      M.N_cov(g,:,p) = round(fracs .* cv);
    end
  end
  fprintf('  Total after:  %d\n', sum(sum(M.N_cov(:,M.TOT,:))));
end

end
