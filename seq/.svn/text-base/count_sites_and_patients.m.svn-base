function M = count_sites_and_patients(M,gta)

  if ~exist('gta','var'), gta = 1:slength(M.gene); end

  demand_fields(M,{'gene','mut'});
  if ~isfield(M,'ng'), M.ng = slength(M.gene); end
  flag1=false;
  if isfield(M.gene,'gene') && ~isfield(M.gene,'name'), M.gene = rename_field(M.gene,'gene','name'); flag1=true; end

  fprintf('Finding number of unique sites, unique patients: ');

  M.gene.npat = zeros(M.ng,1);
  M.gene.nsite = zeros(M.ng,1);
  for g=gta, if ~mod(g,1000), fprintf('%d/%d ',g,M.ng); end
    if ~isfield(M.gene,'genes')   % i.e. if not dealing with genesets
      if isnumeric(M.mut.gene)
        if isfield(M,'use_nonsilent')
          idx = M.use_nonsilent(M.mut.gene(M.use_nonsilent)==g);
        else
          idx = find(M.mut.gene==g);
        end
      else
        if ~isfield(M.mut,'gene_idx')
          M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);
        end
        if isfield(M,'use_nonsilent')
          idx = M.use_nonsilent(M.mut.gene_idx(M.use_nonsilent)==g);
        elseif isfield(M.mut,'is_coding') && isfield(M.mut,'is_silent')
          M.mut = make_numeric(M.mut,{'is_coding','is_silent'});
          idx = find(M.mut.gene_idx==g & M.mut.is_coding & ~M.mut.is_silent);
        else
          error('Need help knowing which mutations to count!  Please supply use_nonsilent or is_coding+is_silent');
        end
      end
    else    % if dealing with genesets
      idx = M.gene.mutnos{g};   %% (because a mutation can belong to more than one geneset)
    end
    if ~isempty(idx)
      M.gene.npat(g) = length(unique(M.mut.patient(idx)));
      M.gene.nsite(g) = length(unique(M.mut.start(idx)));
    end
  end, fprintf('\n');

 if flag1, M.gene = rename_field(M.gene,'name','gene'); end
