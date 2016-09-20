function draw_rate_vs_expr_plot(R,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'fontsize',12);
P = impose_default_value(P,'error_bars','stdev');

R = reorder_struct(R,strcmp(R.conservation,'total'));

ncols = find(~strcmp(R.base_context,R.base_context{1}),1)-1;
[u ui uj] = unique(R.base_context);
types = R.base_context(sort(ui));
[u ui uj] = unique(R.expression_level);
exprtypes = R.expression_level(sort(ui));

rate = reshape(R.mutrate,ncols,slength(R)/ncols);
if strcmpi(P.error_bars,'stdev')
  stdev = reshape(R.stdev,ncols,slength(R)/ncols);
  error_high = rate+stdev;
  error_low = rate-stdev;
elseif strcmpi(P.error_bars,'ci')
  error_high = reshape(R.ci_high,ncols,slength(R)/ncols);
  error_low = reshape(R.ci_low,ncols,slength(R)/ncols);
else
  error('uknown P.error_bars');
end

close,figure(1);clf
barweb(rate',error_high',error_low',0.8,types);
ylabel('mutations per million sites','fontsize',P.fontsize);
legend(exprtypes,'location','northeast','interpreter','none');
set(gca,'fontsize',P.fontsize);

