function M = process_3N_breakdown_stats(M)

fprintf('Processing 3N breakdown stats\n');

M.breakdown.frac = -ones(M.ng,24);
for g=1:M.ng
  gi = find(M.breakdown.gene==g & M.breakdown.stats(:,1) ~= -1);
  if ~isempty(gi)
    stats = sum(M.breakdown.stats(gi,:),1);
    tot = sum(stats);
    if tot>0, M.breakdown.frac(g,:) = stats/tot; end
  end
end

% substitute average values for genes without known transcript models

missing = find(M.breakdown.frac(:,1)<0);
not_missing = setdiff(1:M.ng,missing);
avg = mean(M.breakdown.frac(not_missing,:),1);
M.breakdown.frac(missing,:)=repmat(avg,length(missing),1);

% compute average ratios

M.breakdown.sil = sum(M.breakdown.frac(:,1:8),2);     % mean = 0.2295
M.breakdown.mis = sum(M.breakdown.frac(:,9:16),2);    % mean = 0.7304
M.breakdown.non = sum(M.breakdown.frac(:,17:24),2);   % mean = 0.041

end
