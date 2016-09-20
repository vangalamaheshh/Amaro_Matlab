function do_cov1(samples)

if ~iscell(samples), samples = {samples}; end

get_global_coverage_stats(samples,'global_coverage_by_zone.txt','zone');
