function fh_InsertSizeStatsReport(individual_name,tumor_isz,normal_isz,tumor_lanelist,normal_lanelist,lane_blacklist)
% fh_InsertSizeStatsReport(individual_name,tumor_isz,normal_isz,tumor_lanelist,normal_lanelist,lane_blacklist)
%
% writes the following files to the cwd:
%   *_insert_size_distributions.png
%   *_insert_size_by_flowcell.png
%   *_insert_size_scatterplot.png
%   *_insert_size_QC_report.txt        * = <individual_name>
%   report.html
%   num.mixups.txt   = single line telling count of possible mixups
%
% Mike Lawrence 2009

if nargin<5
  error('input parameters required: individual_name, tumor_isz, normal_isz, tumor_lanelist, normal_lanelist, [lane_blacklist]');
end

% load and process data
P=[];
if exist('lane_blacklist','var'), P.lane_blacklist = lane_blacklist; end
fprintf('Loading TUMOR: '); T = load_isz_file(tumor_isz, tumor_lanelist, P);
fprintf('Loading NORMAL: '); N = load_isz_file(normal_isz, normal_lanelist, P);

fprintf('Generating plots\n');

% raw distributions
plot_TN_isz_distribs(T,N,individual_name);
print('-dpng','-r100',[individual_name '_insert_size_distributions.png']);

% by flowcell
plot_insert_size_by_flowcell(T,N,[],nan,[],individual_name)
print('-dpng','-r100',[individual_name '_insert_size_by_flowcell.png']);

% scatterplot
scatter_width_vs_mean(T,N,individual_name)
print('-dpng','-r100',[individual_name '_insert_size_scatterplot.png']);

% QC report: textfile
R = generate_isz_QC_report(T,N,individual_name);
save_struct(R,[individual_name '_insert_size_QC_report.txt'],[],true);

% QC report: HTML
H = convert_isz_QC_report_to_html(R,individual_name);
save_textfile(H,'report.html');

% n_mixups.txt
n_mixups = sum(R.is_mixup & ~R.is_blacklisted);
save_textfile(sprintf('%d\n',n_mixups),'num.mixups.txt');

% done
fprintf('Report complete.\n');
close all;
