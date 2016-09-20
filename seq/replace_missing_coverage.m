function M = replace_missing_coverage(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'coverage_replace_cutoff',0);
P=impose_default_value(P,'impose_average_cutoff',false);

fprintf('Replacing missing coverage\n');
fprintf('  Total before: %d\n', sum(sum(M.N_cov(:,M.TOT,:))));

Ntot_cov = sum(M.N_cov,3);
f_cov = Ntot_cov ./ (M.np * M.N_terr);
ftot_cov = f_cov(:,M.TOT);
avg_ftot_cov = mean(ftot_cov(ftot_cov>0));        % 0.74
if P.impose_average_cutoff
  avg_ftot_cov = P.impose_average_cutoff;
end

fprintf('  Using average gene coverage of %0.2f\n',avg_ftot_cov);

low_cov = find(ftot_cov <= P.coverage_replace_cutoff);

fprintf('  Number of genes below replacement cutoff: %d\n',length(low_cov));

N_est = round(M.N_terr * avg_ftot_cov);
N_est(:,M.TOT) = sum(N_est(:,1:M.TOT-1-M.NUM_INDEL_CLASSES),2);        % otherwise TOT has rounding errors +-3

% replace coverage

M.N_cov(low_cov,:,:) = repmat(N_est(low_cov,:),[1 1 M.np]);
fprintf('  Total after:  %d\n', sum(sum(M.N_cov(:,M.TOT,:))));

