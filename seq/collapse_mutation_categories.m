function M = collapse_mutation_categories(M,P)
% OBSOLETE

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'collapse_C_and_G',true);
P=impose_default_value(P,'collapse_A_and_T',true);
P=impose_default_value(P,'separate_inframes_and_frameshifts',false);

fprintf('Collapsing mutation categories\n');

  n_missense_tmp = M.n_missense;
  n_silent_tmp = M.n_silent;
  n_nonsense_tmp = M.n_nonsense;
  n_splice_tmp = M.n_splice;
  n_indel_tmp = M.n_indel;
  n_nonsilent_tmp = M.n_nonsilent;
  N_cov_tmp = M.N_cov;
  N_terr_tmp = M.N_terr;

  M.n_missense = [];
  M.n_silent = [];
  M.n_nonsense = [];
  M.n_splice = [];
  M.n_indel = [];
  M.n_nonsilent = [];
  M.N_cov = [];
  M.N_terr = [];

  if P.collapse_C_and_G
    M.mutclass = { ...
      'CpG' ...
      'C+G' ...
    };
    agg = { [3 6] [4 5 7 8] };
  else
    M.mutclass = { ...
      'CpG' ...
      'C' ...
      'G' ...
    };
    agg = { [3 6] [4 5] [7 8] };
  end

  if P.collapse_A_and_T
    M.mutclass = { M.mutclass{:} ...
      'A+T' ...
    };
    agg = { agg{:} [1 2] };
  else
    M.mutclass = { M.mutclass{:} ...
      'A' ...
      'T' ...
    };
    agg = { agg{:} [1] [2] };
  end

  if P.separate_inframes_and_frameshifts
    M.mutclass = { M.mutclass{:} ...
      'InfID' ...
      'Fshft' ...
      'Total' ...
    };
    M.NUM_INDEL_CLASSES = 2;
    agg = { agg{:} [9 11] [10 12 13] M.TOT};
  else
    M.mutclass = { M.mutclass{:} ...
      'Indel' ...
      'Total' ...
     };
    M.NUM_INDEL_CLASSES = 1;
    agg = { agg{:} [9 10 11 12 13] M.TOT};
  end

  M.TOT = length(M.mutclass);

  for a=1:length(agg)
    M.n_missense(:,a,:) = sum(n_missense_tmp(:,agg{a},:),2);
    M.n_silent(:,a,:) = sum(n_silent_tmp(:,agg{a},:),2);
    M.n_nonsense(:,a,:) = sum(n_nonsense_tmp(:,agg{a},:),2);
    M.n_splice(:,a,:) = sum(n_splice_tmp(:,agg{a},:),2);
    M.n_indel(:,a,:) = sum(n_indel_tmp(:,agg{a},:),2);
    M.n_nonsilent(:,a,:) = sum(n_nonsilent_tmp(:,agg{a},:),2);
    M.N_cov(:,a,:) = sum(N_cov_tmp(:,agg{a},:),2);
    M.N_terr(:,a,:) = sum(N_terr_tmp(:,agg{a},:),2);
  end

  % restore correct number for indel bases at risk after summing

  M = define_indel_bases_at_risk(M);

  % collapse breakdown.frac

  nc = length(agg)-M.NUM_INDEL_CLASSES-1;
  bftmp = M.breakdown.frac;
  M.breakdown.frac = zeros(M.ng,3*nc);
  for a=1:nc
    M.breakdown.frac(:,a) = sum(bftmp(:,agg{a}),2);            %sil
    M.breakdown.frac(:,nc+a) = sum(bftmp(:,8+agg{a}),2);       %mis
    M.breakdown.frac(:,2*nc+a) = sum(bftmp(:,16+agg{a}),2);    %non
  end

end
