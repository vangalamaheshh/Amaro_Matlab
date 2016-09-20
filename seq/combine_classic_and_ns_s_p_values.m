function M = combine_classic_and_ns_s_p_values(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'sigmoid_parameter_alpha',0.5);
P=impose_default_value(P,'sigmoid_parameter_beta',6);
P=impose_default_value(P,'pval_cutoff',1e-11);

demand_fields(M,'gene');
demand_fields(M.gene,{'pval_classic','pval_ns_s'});

if ~isnumeric(M.gene.pval_classic) || ~isnumeric(M.gene.pval_ns_s)
  error('~isnumeric(M.gene.pval_classic) || ~isnumeric(M.gene.pval_ns_s)');
end

alpha = P.sigmoid_parameter_alpha;
beta = P.sigmoid_parameter_beta;

M.gene.pval = nan(M.ng,1);
for i=1:M.ng
%  tmp = (M.gene.pval_ns_s(i)/alpha)^beta;
%  f = 1-(tmp/(tmp+1));
  f = 1 / (1 + ((M.gene.pval_ns_s(i)/alpha)^beta));
  M.gene.pval(i) = (f*M.gene.pval_classic(i)) + ((1-f)*M.gene.pval_ns_s(i));
end

% make sure no p>1 !
M.gene.pval = min(1,M.gene.pval);

% fix negative and unreliable very-small values
M.gene.pval_lessthan_flag = (M.gene.pval <= P.pval_cutoff);
M.gene.pval(M.gene.pval_lessthan_flag) = P.pval_cutoff;

