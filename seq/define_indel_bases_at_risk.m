function M = define_indel_bases_at_risk(M)

% every base is at risk for an indel

M.N_cov(:,M.TOT-M.NUM_INDEL_CLASSES:M.TOT-1,:) = repmat(M.N_cov(:,M.TOT,:),1,M.NUM_INDEL_CLASSES);

end
