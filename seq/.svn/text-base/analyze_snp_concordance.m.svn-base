function [r] = analyze_snp_concordance(snpconc, lane_blacklist, output_prefix)

a = load_struct(snpconc, '%s%s%d%s%s%s%s%s%f%d%d%d%d%d%d');

[u,ui,uj]=unique(a.lane);                                        
qq=accumcols(uj,[a.major_allele_counts a.minor_allele_counts]);

r.lane = u;
r.major_allele_counts = qq(:,1);
r.minor_allele_counts = qq(:,2);
r.concordance = r.minor_allele_counts./(r.major_allele_counts + r.minor_allele_counts);

% identify blacklisted lanes
if strcmpi(lane_blacklist,'none')
  BL = [];
else
  BL = load_lines(lane_blacklist);
end
r.is_blacklisted = is_blacklisted(r.lane,BL);

save_struct(r, [output_prefix '.snp_concordance_summary.txt']);
