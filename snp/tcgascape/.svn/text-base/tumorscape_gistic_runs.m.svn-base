function tumorscape_gistic_runs(input_parameter_file,run_dir)
% TUMORSCAPE_GISTIC_RUNS launch GISTIC runs for tumorscape

% do gistic runs for tumorscape
% 16-Apr-2012 use "old method" with a different segfile for each disease
% this is a run on hg19 data normalized with tangent from the "2.1"
% pipeline. This run does not write gene files and uses array list from
% TCGA/GDAC firehose run on 21-Mar-2012 to filter samples

if ~exist('run_dir','var') || isempty(run_dir)
    run_dir = [pwd filesep];
end

% read tumorscape input parameter file (XML document)
TSP = read_tumorscape_params(input_parameter_file);

% method for determining disease type
%   true  = use a sample info on single pan-cancer seg file
%   false = use separate segfiles for each disease, mapped by
%           tcga_segfgiles() !!!TODO figure out a scheme for this
sample_info_method = true;

output_dir = run_dir; %!'/xchip/gistic/tcgascape/tcgascape_120416/';

% cache file for cleaned D
cached_clean_master = [output_dir,'Dmaster.mat'];

% create output directory
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%{
% load input segmented copy number file mapping
segfile_map = tcgascape_segfiles();
%}

if exist(cached_clean_master,'file')
    % load cleaned data from cache file
    D = load_D(cached_clean_master);
else
    % create D from segments
    opts = struct('use_segarray',true,...
                  'compress_markers',true,...
                   'chr_field',false,...
                   'marker_field',false );
               
    if sample_info_method
        % use sample info file to map TCGA id to gcmtype
        SI = read_sample_info_file(TSP.files.sampleinfo);
        gcm_map = containers.Map({SI.sample_name},lower({SI.tcga_experiment}));
        D = make_D_from_seg(TSP.files.segdata,TSP.files.markers,opts);
        %%%
        D = reorder_D_cols(D,isKey(gcm_map,D.sdesc));
        %%%
        D.sis = struct('gcmtype',cellfun(@(x) gcm_map(x),D.sdesc,'UniformOutput',false)); 
    else
        assert(false); %!!! this no longer works
        % load all the data into a master D struct,
        % disease inferred from input file
        D = load_D_tree('all_cancers',TSP.cancer_treegen,TSP.segfile_map,...
                        TSP.files.markers,opts);
    end
    % if array list provided, use it to filter samples
    if ~isempty(TSP.files.arraylist)
        if exist(TSP.files.arraylist,'file')
            AL = read_array_list_file(TSP.files.arraylist);
            [~,~,dx] = match_string_sets_hash({AL.array},D.sdesc);
            D = reorder_D_cols(D,dx);
        else
            error('array list file not found'); %!!! improve diagnostic
        end
    end
    % remove X/Y chromosomes, markers with NaNs, and smooth
    D = clean_gistic_input(D,TSP.files.cnvlist,TSP.gistic.remove_X,...
                           TSP.gistic.join_segment_size);
    % save cached copy of all data
    save_D(cached_clean_master,D,'-v7.3');
end

%% save D structs and segfiles for IGV
gistic_dir = [output_dir 'gistic_inputs/'];
gistic_results_dir = [output_dir 'gistic_results/'];
segfile_dir = [output_dir 'igvfiles/'];

if ~exist(gistic_dir,'dir');mkdir(gistic_dir);end
if ~exist(gistic_results_dir,'dir');mkdir(gistic_results_dir);end
if ~exist(segfile_dir,'dir');mkdir(segfile_dir);end

splitsave_D(D,'all_cancers',TSP.cancer_treegen,gistic_dir,segfile_dir,...
            TSP.cancer_namemap);

%% run gistic on all disease categories
run_gistic_on_all(gistic_dir,gistic_results_dir,TSP.files.refgene,...
                  TSP.gistic,TSP.gistic_cpus,TSP.min_samples);





