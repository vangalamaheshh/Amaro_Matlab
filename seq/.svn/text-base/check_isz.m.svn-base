function check_isz(samples,P,T,N)
% check_isz(samples,P)
%
% Mike Lawrence 2009-09-10

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'InsertSizeByLane_output_dir_suffix','isz');
P = impose_default_value(P,'InsertSizeByLane_output_file_extension','isz');
P = impose_default_value(P,'InsertSizeByLane_java_classname','InsertSizeByLane');
P = impose_default_value(P,'InsertSizeByLane_code','ISZ');
P = impose_default_value(P,'InsertSizeByLane_max_isz','2000');

if ~exist('T','var') | ~exist('N','var')
  [T N] = load_isz_data(samples,P);
  if length(samples)==1
    tmp = T; T = {tmp};
    tmp = N; N = {tmp};
  end
end

for i=1:length(samples)

  fprintf('\nSAMPLE %d/%d = %s\n',i,length(samples),samples{i});

  % distributions
  plot_TN_isz_distribs(T{i},N{i},samples{i});
  fprintf('  Plot 1/3.  Type "return" for next plot.\n');
  keyboard

  % insert size by flowcell
  plot_insert_size_by_flowcell(T{i},N{i},[],nan,[],samples{i})
  fprintf('  Plot 2/3.  Type "return" for next plot.\n');
  keyboard

  % peak width vs. adjusted mean
  scatter_width_vs_mean(T{i},N{i},samples{i})
  fprintf('  Plot 3/3.  Type "return" for next sample.\n');
  keyboard
end

fprintf('\nType "return" to end.\n');
keyboard
