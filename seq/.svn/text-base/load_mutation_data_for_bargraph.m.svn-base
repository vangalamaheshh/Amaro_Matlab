function D = load_mutation_data_for_bargraph(M,P)
% OBSOLETE: now done within generate_coverage_plot_and_bargraphs

if ~exist('P','var'), P=[]; end

D = [];

if isfield(M,'maffile') && isfield(M,'covfile') && isfield(M,'patfile')

  demand_file(M.maffile);
  demand_file(M.covfile);
  demand_file(M.patfile);

  % input from list of files
  I=M;
  M=[];
  M.patient = load_struct(I.patfile);
  if isfield(I,'ordfile')
    tmp = load(I.ordfile);
    if length(tmp.snameidx)~=slength(M.patient), error('ordfile and patfile different lengths'); end
    if isnumeric(tmp.snameidx)
      ord = tmp.snameidx;
    else
      ord = listmap(tmp.snameidx,M.patient.dir);
    end
    n = sum(isnan(ord));
    if n>0, error('%d patfile patients not in ordfile',n); end
    M.patient = reorder_struct(M.patient,ord);
  end
  tmp = load(I.covfile);
  if isfield(tmp,'COV')
    if ~isfield(tmp.COV,'sample') && isfield(tmp.COV,'patient')
      tmp.COV = rename_field(tmp.COV,'patient','sample');
    end
    if isfield(M.patient,'dir')
      idx = listmap(M.patient.dir,tmp.COV.sample.dir);
    else
      idx = listmap(M.patient.name,tmp.COV.sample.name);
    end
    n = sum(isnan(idx));
    if n>0, error('%d patfile patients not in covfile',n); end
    N_tot = sum(sum(tmp.COV.cov(:,idx,:),3),1)';
  elseif isfield(tmp,'C1')
    if isfield(M.patient,'dir')
      idx = listmap(M.patient.dir,tmp.C1.sample.dir);
    else
      idx = listmap(M.patient.name,tmp.C1.sample.name);
    end
    D = reorder_struct(tmp.C1.sample,idx);
    n = sum(isnan(idx));
    if n>0, error('%d patfile patients not in covfile',n); end
    if isfield(tmp.C1,'totcov')
      N_tot = sum(tmp.C1.totcov(:,idx))';
    elseif isfield(tmp.C1,'cov')
      N_tot = sum(sum(tmp.C1.cov(:,idx,:),3),1)';
    else
      error('C1 contains neither cov nor totcov');
    end
  else
    error('covfile contains neither COV nor C1');
  end
  tmp = load_struct(I.maffile);

  sil = grep('Synonymous|Silent',tmp.type,1);
  non = grep('Mis|Nons|Read-through|Splice',tmp.type,1);
  if isfield(M.patient,'dir')
    sno = listmap(tmp.sample_dir,M.patient.dir);
  else
    sno = listmap(tmp.patient,M.patient.name);
  end
  nsil_tot = histc(sno(sil),1:slength(M.patient));
  nnon_tot = histc(sno(non),1:slength(M.patient));

elseif isfield(M,'n_silent') && isfield(M,'n_nonsilent') && isfield(M,'N_cov')
  % input from TCGA-style data structure
  P=impose_default_value(P,'sort_order',1:M.np);
  ord = P.sort_order;
  nsil_tot = sum(squeeze(M.n_silent(:,M.TOT,ord)))';
  nnon_tot = sum(squeeze(M.n_nonsilent(:,M.TOT,ord)))';
  N_tot = sum(squeeze(M.N_cov(:,M.TOT,ord)))';
  M.patient = reorder_struct(M.patient,ord);
else
  error('Unknown input format');
end

if isfield(M.patient,'short')
  names = M.patient.short;
else
  names = M.patient.name;
  names = regexprep(names,'.*-(\d+)(-.*)?','$1');
end

rate_sil = nsil_tot ./ N_tot;
rate_non = nnon_tot ./ N_tot;

D.display_name = names;
D.nnon_tot = nnon_tot;
D.nsil_tot = nsil_tot;
D.N_tot = N_tot;
D.rate_sil = rate_sil;
D.rate_non = rate_non;
