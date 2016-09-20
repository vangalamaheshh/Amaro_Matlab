function report_insert_size_stats(sample,P)
% report_insert_size_stats(sample,P)
%
% Mike Lawrence 2009-07-29

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Does not support multiple samples at this time.'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'InsertSizeByLane_output_dir_suffix','isz');
P = impose_default_value(P,'insert_size_distributions_filename','insert_size_distributions.pdf');
P = impose_default_value(P,'insert_size_by_flowcell_filename','insert_size_by_flowcell.pdf');
P = impose_default_value(P,'insert_size_scatterplot_filename','insert_size_scatterplot.pdf');

basedir = '/xchip/tcga_scratch/lawrence';

fprintf('Loading insert-size statistics\n');
[T N] = load_TN_perlane_data(sample,P);

fprintf('Generating plots\n');

% raw distributions
plot_TN_isz_distribs(T,N);
fname = [basedir '/' sample '/' P.insert_size_distributions_filename];
print('-dpdf',fname);

% by flowcell
plot_insert_size_by_flowcell(T,N,[],nan)
fname = [basedir '/' sample '/' P.insert_size_by_flowcell_filename];
print('-dpdf',fname);

% scatterplot
scatter_width_vs_mean(T,N)
fname = [basedir '/' sample '/' P.insert_size_scatterplot_filename];
print('-dpdf',fname);

% done
close all;
