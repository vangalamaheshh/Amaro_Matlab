function [D] = get_filtered_dbSNP_stats_v2(workspace, individual_set, varargin)
%% Creates a table describing the fraction of filtered dbSNPs in 
%% each sample 
MAX_NUM_ARGS = 4;
if nargin > MAX_NUM_ARGS, error('Too many arguments! A max of four is allowed.\n'), end

flag = varargin{1};
if nargin == MAX_NUM_ARGS
  out_file = varargin{2};
end

subpath = fullfile('/xchip', 'cga1', 'firehose_output', workspace, 'Individual_Set', individual_set);
f = direc(subpath);
patients = regexprep(f, strcat(subpath, '/'), '');
if strcmp(flag, 'filter')
  patients_suffix = cellfun(@strcat, patients, repmat({'.normal_filtered.call_stats.txt'}, length(patients), 1), 'UniformOutput', false);
elseif strcmp(flag, 'callstats')
  patients_suffix = cellfun(@strcat, patients, repmat({'.call_stats.txt'}, length(patients), 1), 'UniformOutput', false);
else 
  error('Either ''filter'' or ''callstats'' are allowed for the FH output subdirectory path retrieval.\n');
end
%% Find capture paths

if strcmp(flag, 'filter');
  paths_to_annotated = cellfun(@fullfile, f, repmat({'capture'}, length(f), 1), ...
                               repmat({'mut'}, length(f), 1),  repmat({'normal_filtered'}, length(f), 1), patients_suffix, ...
                               'UniformOutput', false);
else
  paths_to_annotated = cellfun(@fullfile, f, repmat({'capture'}, length(f), 1), ...
                               repmat({'mut'}, length(f), 1),  repmat({'calls'}, length(f), 1), patients_suffix, ...
                               'UniformOutput', false);
end
idx = cellfun(@file_exist, paths_to_annotated);
paths_to_annotated_capture = paths_to_annotated(find(idx));
patients_capture = patients(find(idx));


%% Find WGS paths 
if strcmp(flag, 'filter');
  paths_to_annotated = cellfun(@fullfile, f, repmat({'wgs'}, length(f), 1), ...
                               repmat({'mut'}, length(f), 1),  repmat({'normal_filtered'}, length(f), 1), patients_suffix, ...
                               'UniformOutput', false);
else
  paths_to_annotated = cellfun(@fullfile, f, repmat({'wgs'}, length(f), 1), ...
                               repmat({'mut'}, length(f), 1),  repmat({'calls'}, length(f), 1), patients_suffix, ...
                               'UniformOutput', false);
end
  

idx = cellfun(@file_exist, paths_to_annotated);
paths_to_annotated_wgs = paths_to_annotated(find(idx));
patients_wgs = patients(find(idx));

patients = [patients_capture, patients_wgs];
paths_to_annotated = [paths_to_annotated_capture, paths_to_annotated_wgs];

%keyboard
%% Fill out output struct 

D = []; 
D.individual_id = patients; 
D.total = nan(slength(D), 1);
D.DBSNP = nan(slength(D), 1);
D.fraction_DBSNP = nan(slength(D), 1);
D.percentile_DBSNP_AF = nan(slength(D), 1);
D.num_dbsnp_seen_in_other_normal = nan(slength(D), 1);
for i = 1:slength(D)
  disp(i)
  try 
    maf = load_struct(paths_to_annotated{i});
  catch
    disp(skipping)
    disp(D.individual_id)
    continue 
  end
  try
    if strcmp(flag, 'filter')
      idx = find(strcmp(maf.dbsnp_site, 'DBSNP') ...
                 & strcmp(maf.failure_reasons, 'seen_in_normal_panel') ... 
                 & strcmp(maf.normal_f, '0')); 
    else 
      idx = find(strcmp(maf.dbsnp_site, 'DBSNP') ...
                         & strcmp(maf.failure_reasons, 'seen_in_panel_of_normals') ...        
                         & strcmp(maf.normal_f, '0'));
    end
  catch
    keyboard
    continue
  end
  if ~isempty(idx) 
    D.num_dbsnp_seen_in_other_normal(i) = length(idx); 
    maf = make_numeric(maf, 'tumor_f');
    p = prctile(maf.tumor_f(idx), 70); 
%    hist(maf.tumor_f(idx), 100); 
%    xlim([0 1]);
%    print_to_file(['/xchip/cga1/petar/MM/data_qc/perc_95/', D.individual_id{i}, '.dbsnp_hist.png']);
    D.DBSNP(i) = length(idx);
    D.fraction_DBSNP(i) = length(idx)/slength(maf); 
    D.percentile_DBSNP_AF(i) = p; 
    D.total(i) = slength(maf);
  end 
 % keyboard
  %D.DBSNP(i) = b(idx); 
  %D.nonDBSNP(i) = sum(b) - D.DBSNP(i);
  %D.ratio(i) = D.DBSNP(i)./sum(b);
end 

if nargin == 4 
  save_struct(D, out_file);
end 



 %   normals = cellfun(@splitstr, maf.nomal_samples_called_names(idx), repmat({','}, length(idx), 1), ...
  %                  'UniformOutput', false);
%    normal_paths = cell(length(normals), 1);
%    for j = 1:length(normals) 
%      normal_paths{j} = cellfun(@regexprep, repmat({paths_to_annotated{i}}, length(normals{j}), 1), ...
%                        repmat({[D.individual_id{i}, '.normal_filtered.call_stats.txt']}, length(normals{j}), 1), normals{j}, 'UniformOutput', false); 
%      normal_paths{j}  = regexprep(normal_paths{j}, '.cleaned', '');
%      normal_paths{j}  = regexprep(normal_paths{j}, '-Normal','-Normal.cleaned.call_stats.txt');
%      idx = find(demand_file(normal_paths{j}));
%      N = load_structs(normal_paths{j}(idx)); 
 
%    end
