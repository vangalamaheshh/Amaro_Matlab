function [G m v] = analyze_mutrate_covariates_and_quick_mutsig_test(varargin)

if nargin>1 && isstruct(varargin{end}) && ~isfield(varargin{end},'val')
  vargs = varargin(1:end-1);
  P = varargin{end};
else
  vargs = varargin;
  P = [];
end

G = vargs{1};
ng = slength(G);

P = impose_default_value(P,'genes_to_analyze',1:ng);
P = impose_default_value(P,'report_top_n_mutsig_genes',min(40,length(P.genes_to_analyze)));
P = impose_default_value(P,'generate_plot',false);

tmp = P.generate_plot;
P.generate_plot = false;  % (wait til after adding MutSig ranks)

if nargout==3
  [G m v] = analyze_mutrate_covariates(vargs{:},P);
else
  G = analyze_mutrate_covariates(vargs{:},P);
end

G = quick_mutsig_test(G,P);

if tmp
  plot_Ffit_vs_Fobs(G,P);
end
