function plot_allelic_fracs_with_CN(workspace, individual_set, out_dir, varargin)
% Plots allelic fractions with beta distributions and if CN data is available,
% with colored histograms showing the relative CN log ratio of each mutation on the histogram. 
% If SIF file is unavailable, will only do allelic fraction graphs 
% Parameters: FH workspace and individual set, and optionally a SIF file linking 
% each FH individual to the respective segmentation file. 
%
%
%

%% Retrieve ready .maf files from FH output directory 
subpath = fullfile('/xchip', 'cga1', 'firehose_output', workspace, 'Individual_Set', individual_set);
f = direc(subpath);
patients = regexprep(f, strcat(subpath, '/'), '');
patients_suffix = cellfun(@strcat, patients, repmat({'-Tumor.maf.annotated'}, length(patients), 1), 'UniformOutput', false);
paths_to_annotated = cellfun(@fullfile, f, repmat({'capture'}, length(f), 1), ...
                             repmat({'mut'}, length(f), 1),  repmat({'annotated'}, length(f), 1), patients_suffix, ...
                             'UniformOutput', false);
idx = cellfun(@file_exist, paths_to_annotated);
paths_to_annotated = paths_to_annotated(find(idx));
patients = patients(find(idx));

%% Process rest of parameters 
ensure_dir_exists(out_dir);
outstem = out_dir;
if ~isempty(varargin)
  seg_path = varargin{1};
  S = load_struct_noheader(seg_path);
  S = reorder_struct(S, ~strcmp(S.col2, ''));
  S.Individual_Id= S.col1;
  S.snp_array = S.col2;
else 
  S = []; 
  S.Individual_Id = cell(1,1); 
  S.snp_array = cell(1,1); 
end

f = paths_to_annotated;
%% Create plots and save to output directory 

for i = 1:length(f)
    path = f{i};
    a = load_struct(path);
    if isempty(a.Hugo_Symbol), continue, end
    disp(patients{i});
    if isfield(a, 'tumor_f')
      a = parse_in(a, {'tumor_f', 't_ref_count', 't_alt_count'}, '.*', {'i_tumor_f', 'i_t_ref_count', 'i_t_alt_count'});
    end
    a = make_numeric(a, {'i_tumor_f', 'i_t_ref_count', 'i_t_alt_count', 'Chromosome', 'Start_position'});
    max_y = max(histc(a.i_tumor_f, 0:0.01:1)) + 1;
%    keyboard
    if ~ismember(patients{i}, S.Individual_Id)
      betas = nan(slength(a), length(0:0.01:1));
      for j = 1:slength(a)
        x = a.i_t_alt_count(j);
        n = a.i_t_ref_count(j) + a.i_t_alt_count(j);
        beta = betapdf(0:0.01:1, x+1,n-x+1);
        betas(j,:) = beta;
      end

      subplot(3, 1, 1);
      title(patients{i});
      plot(0:0.01:1,betas)
      xlim([0 1]);
      subplot(3, 1, 2);
      plot(0:0.01:1,sum(betas,1)./slength(a));
      xlim([0 1]);
      subplot(3, 1, 3);
      hist(a.i_tumor_f, 50);
      xlim([0 1]);            
      ylim([0 max_y+1]);
      print_to_file([outstem '/' patients{i} '.all_figs.png']);
      continue;
    end

    idx = find(strcmp(patients{i}, S.Individual_Id));
    array = S.snp_array(idx);
    plate = splitstr(char(array), '_');
    plate = plate{1};
    path = ['/xchip/cga2/analysis_pipeline_plates2/cn/latest_runs/production_gap/' ...
        plate '/BySample/seg/' char(array) '.seg.data.txt'];
    seg = load_struct(path);
    seg = make_numeric(seg, {'chrom', 'locstart', 'locend', 'segmean'});
    %  keyboard
    betas = nan(slength(a), length(0:0.01:1));
    for j = 1:slength(a)
        x = a.i_t_alt_count(j);
        n = a.i_t_ref_count(j) + a.i_t_alt_count(j);
        beta = betapdf(0:0.01:1, x+1,n-x+1);
        betas(j,:) = beta;
    end
    color_in = nan(slength(a), 4);
    color_in(:, 3:4) = repmat(0, slength(a), 2);
    seg_muts = nan(slength(a), 1);
    for j = 1:slength(a)
        idx = find(a.Chromosome(j) == seg.chrom & a.Start_position(j) >= seg.locstart & a.Start_position(j) <= seg.locend);
        if isempty(idx) continue; end
        segmean = seg.segmean(idx);
        seg_muts(j) = segmean;
        color_in(j, 1) = a.i_tumor_f(j);
        if segmean > -0.2 & segmean < 0.2
            color_in(j, 2:4) = [1 1 1];
        elseif segmean <= -0.2
            color_in(j, 2:4) = [0 0 1];
        elseif segmean >= 0.2
            color_in(j, 2:4) = [1 0 0];
        end
    end
    idx = find(~isnan(color_in(:, 1)));
    color_in = color_in(idx, :);
    seg_muts = seg_muts(idx);

    subplot(5, 1, 1);
    title(patients{i});
    plot(0:0.01:1,betas)
    subplot(5, 1, 2);
    plot(0:0.01:1,sum(betas,1)./slength(a));
    subplot(5, 1, 3);
    hist(a.i_tumor_f, 100);
    xlim([0 1]);
    ylim([0 max_y+1]);
%    keyboard
    subplot(5, 1, 4);
    colored_hist(color_in(:, 1), color_in(:, 2:4), 100);
    xlim([0 1]);
    ylim([0 max_y+1]);
    subplot(5, 1, 5);
    hist(seg_muts, 50);
    xlim([-2 2]);
    print_to_file([outstem '/' patients{i} '.all_figs.png']);
end
%     catc









