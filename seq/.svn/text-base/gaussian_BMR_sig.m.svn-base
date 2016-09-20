function P_sig = gaussian_BMR_sig(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'permute_count',1000);
P=impose_default_value(P,'Q_threshold',0.1);
P=impose_default_value(P,'gaussian_method',1);
P=impose_default_value(P,'gaussian_center',M.mutrate.tot.hat);
P=impose_default_value(P,'gaussian_std',M.mutrate.tot.stdev);

BMR = random('normal',P.gaussian_center,P.gaussian_std,P.permute_count,1);
P.mutation_rate_to_use = 'manual';

if P.gaussian_method == 1

  sig_count = zeros(M.ng,1);

  for run_no = 1:P.permute_count
    fprintf('Run %d of %d\n', run_no, P.permute_count);
    P.manual_mutation_rate = BMR(run_no);  
    M = find_significant_genes(M,P);
    sig_genes = find(M.Q < P.Q_threshold);
    sig_count(sig_genes) = sig_count(sig_genes) + 1;
  end

  P_sig = sig_count / P.permute_count;


elseif P.gaussian_method == 2

  P_tot = zeros(M.ng,1);

  for run_no = 1:P.permute_count
    fprintf('Run %d of %d\n', run_no, P.permute_count);
    P.manual_mutation_rate = BMR(run_no);
    M = find_significant_genes(M,P);
    P_tot = P_tot + M.Prob;
  end

  P_avg = P_tot / P.permute_count;
  P_sig = P_avg;

end

