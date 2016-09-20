function R = analyze_gene_gene_correlations(M,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'mutsig_significance_threshold',0.1);
P = impose_default_value(P,'max_genes_to_consider',50);
P = impose_default_value(P,'gene_gene_correlations_num_permutations',10000);
P = impose_default_value(P,'gene_gene_correlations_report_filename','*required*');
P = impose_default_value(P,'correlation_significance_threshold',0.05);

if ~isfield(M.gene,'qval')
  error('Requires significance analysis to have been run');
end

R = [];
R.gene1 = {};
R.gene2 = {};
R.relationship = {};
R.p = [];
R.q = [];

sig_genes = find(M.gene.qval<=P.mutsig_significance_threshold);
if length(sig_genes)>P.max_genes_to_consider
  fprintf('Truncating list of significant genes after the first %d genes\n',P.max_genes_to_consider);
  [tmp ord] = sort(M.gene.qval);
  sig_genes = ord(1:P.max_genes_to_consider);
end

if length(sig_genes)<2
  fprintf('Need at least 2 significantly mutated genes to run algorithm.\n');
else

  if isnumeric(M.mut.gene) && isnumeric(M.mut.patient)
    T = [M.mut.gene(M.use_nonsilent) M.mut.patient(M.use_nonsilent)];
  elseif isfield(M.mut,'gene_idx') && isfield(M.mut,'pat_idx')
    T = [M.mut.gene_idx(M.use_nonsilent) M.mut.pat_idx(M.use_nonsilent)];
  else
    error('waah!');
  end

  nruns = P.gene_gene_correlations_num_permutations;

  npairs = (length(sig_genes) * (length(sig_genes)-1)) / 2;
  R.gene1 = cell(npairs,1);
  R.gene2 = cell(npairs,1);
  R.relationship = cell(npairs,1);
  R.p = nan(npairs,1);

  ridx = 1;
  for g1=1:length(sig_genes)-1
    for g2=g1+1:length(sig_genes)
      g1idx = sig_genes(g1);
      g2idx = sig_genes(g2);
      fprintf('%s vs. %s',M.gene.name{g1idx},M.gene.name{g2idx});
      m1_obs = find(T(:,1)==g1idx);
      m2_obs = find(T(:,1)==g2idx);
      pc_obs = intersect(T(m1_obs,2),T(m2_obs,2));
      Nc = sum(ismember(T([m1_obs;m2_obs],2),pc_obs));
      Ne = length([m1_obs;m2_obs])-Nc;
      fprintf('\n\tm1 %d   m2 %d   pc %d    Nc %d    Ne %d\n',...
              length(m1_obs),length(m2_obs),length(pc_obs),Nc,Ne);
      c=0;
      len_pc_obs = length(pc_obs);
      pc_ct = 0;
      for i=1:nruns
        Tp = [T(randperm(size(T,1)),1) T(:,2)];
        m1 = find(Tp(:,1)==g1idx);
        m2 = find(Tp(:,1)==g2idx);
        pc = intersect(Tp(m1,2),Tp(m2,2));
        Xc = sum(ismember(Tp([m1;m2],2),pc));
        Xe = length([m1;m2])-Xc;
        c = c + (Xc>=Nc & Xe>=Ne);
        pc_ct = pc_ct + (length(pc) >= len_pc_obs);
%        fprintf('\t\tPermutation %d\n',i);
%        fprintf('\t\tm1 %d   m2 %d   pc %d    Xc %d    Xe %d   result %d\n',...
%              length(m1),length(m2),length(pc),Xc,Xe,double(Xc>=Nc & Xe>=Ne));
      end
      p = c/nruns;
      pc_frac = pc_ct / nruns;
      fprintf('\tp=%f in %d runs      pc_frac = %d\n',p,nruns,pc_frac);
      R.gene1{ridx} = M.gene.name{g1idx};
      R.gene2{ridx} = M.gene.name{g2idx};
      if p<=P.correlation_significance_threshold
        if pc_frac < 0.5
          R.relationship{ridx} = 'correlated';
        else
          R.relationship{ridx} = 'anti-correlated';
        end
      else
        R.relationship{ridx} = 'none';
      end
      R.p(ridx) = p;
      ridx = ridx + 1;
    end    
  end
  R.q = calc_fdr_value(R.p);
  R = sort_struct(R,{'q','p'});
end


save_struct(R,P.gene_gene_correlations_report_filename);


return


%%%% examples for debugging


% mutex
T=[1 1;1 2;1 3;1 4;1 5;2 5;2 7;2 8;2 9; 2 10;...
    3 1;4 2;5 3;6 4; 7 5; 8 6; 9 7; 10 8; 11 9; 12 10;...
   13 1;14 2;15 3;16 4; 17 5; 18 6; 19 7; 20 8; 21 9; 22 10]
g1idx=1;
g2idx=2;

% corr
T=[1 1;1 3;1 5;1 7;1 9;2 1;2 3;2 5;2 7; 2 9;...
    3 1;4 2;5 3;6 4; 7 5; 8 6; 9 7; 10 8; 11 9; 12 10;...
   13 1;14 2;15 3;16 4; 17 5; 18 6; 19 7; 20 8; 21 9; 22 10]
g1idx=1;
g2idx=2;

