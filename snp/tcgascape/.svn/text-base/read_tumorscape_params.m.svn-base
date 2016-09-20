function params = read_tumorscape_params(fname)
%READ_TCGASCAPE_PARAMS read tumorscape run parameters from XML file
typedoc = xmlread(fname);
rootElement = getDocumentElement(typedoc);
root_tag = getTagName(rootElement);
assert(strcmp(root_tag,'tumorscape-parameters'));

nodelist = getElementsByTagName(rootElement,'cancer-root-type');
assert(getLength(nodelist)==1);
cancerTypeElement = item(nodelist,0);

% extract the cancer types
[cancer_tree,cancer_names] = recurse_type(containers.Map,containers.Map,cancerTypeElement);
%gistic = struct('cancer_tree',cancer_tree,'cancer_names',cancer_names);


%% extract GISTIC parameters
nodelist = getElementsByTagName(rootElement,'gistic-run-parameters');
assert(getLength(nodelist)==1);
gisticParameters=item(nodelist,0);

gistic = struct;
gistic.t_amp = get_param(gisticParameters,'t_amp','double');
gistic.t_del = get_param(gisticParameters,'t_del','double');
gistic.broad_len_cutoff = get_param(gisticParameters,'broad_len_cutoff','double');
gistic.cap = get_param(gisticParameters,'cap','double');
gistic.conf_level = get_param(gisticParameters,'conf_level','double');
gistic.join_segment_size = get_param(gisticParameters,'join_segment_size','double');
gistic.qv_thresh = get_param(gisticParameters,'qv_thresh','double');
gistic.ziggs.max_segs_per_sample = get_param(gisticParameters,'max_segs_per_sample','double');
gistic.do_gene_gistic = get_param(gisticParameters,'do_gene_gistic','boolean');
gistic.remove_X = get_param(gisticParameters,'remove_X','boolean');
gistic.arm_peeloff = get_param(gisticParameters,'arm_peeloff','boolean');
gistic.run_broad_analysis = get_param(gisticParameters,'run_broad_analysis','boolean');
gistic.save_seg_data = get_param(gisticParameters,'save_seg_data','boolean');
gistic.write_gene_files = get_param(gisticParameters,'write_gene_files','boolean');
gistic.conserve_disk_space = get_param(gisticParameters,'conserve_disk_space','boolean');
gistic.save_data_files = get_param(gisticParameters,'save_data_files','boolean');
gistic.use_segarray = get_param(gisticParameters,'use_segarray','boolean');

%% extract GISTIC input files
nodelist = getElementsByTagName(rootElement,'gistic-input-files');
assert(getLength(nodelist)==1);
gisticFiles=item(nodelist,0);

files = struct;
files.segdata = get_param(gisticFiles,'segdata');
files.arraylist = get_param(gisticFiles,'arraylist');
files.markers = get_param(gisticFiles,'markers');
% allow multiple CNV lists, separated by whitespace, ',', ';' or '|'
cnv_lists = get_param(gisticFiles,'cnvlist');
files.cnvlist = regexp(cnv_lists,'[\s;,|]+','split');

files.sampleinfo = get_param(gisticFiles,'sampleinfo');
files.refgene = get_param(gisticFiles,'refgene');

params = struct('gistic',gistic,'files',files,'cancer_treegen',cancer_tree,...
                'cancer_namemap',cancer_names);

%% add tumorscape parameters
params.root_disease = get_param(rootElement,'root_disease');
params.min_samples = get_param(rootElement,'min_samples','double');
params.gistic_cpus = get_param(rootElement,'gistic_cpus','double');
params.ht_amp = get_param(rootElement,'ht_amp','double');
params.ht_del = get_param(rootElement,'ht_del','double');
params.run_id = get_param(rootElement,'run_id');
params.igvfile_basedir = get_param(rootElement,'igvfile_basedir');
params.igvfile_baseurl = get_param(rootElement,'igvfile_baseurl');

%% SUBFUNCTION: read a parameter
function value = get_param(xmlElement,attribute,type)

if ~exist('type','var') || isempty(type)
    type='char';
end
charval = char(getAttribute(xmlElement,attribute));
switch type
    case 'char'
        value = charval;
    case 'double'
        value = str2num(charval);
    case 'boolean'
        if strcmp(charval,'true')
            value = true;
        elseif strcmp(charval,'false')
            value = false;
        else
            value = boolean(str2double(charval));
        end
    otherwise
        error('unsupported attribute type');
end
        

%% SUBFUNCTION: parse the cancer types used in the analysis
function [cancer_tree,cancer_name] = recurse_type(cancer_tree,cancer_name,element)
subtypes = getChildNodes(element);
id = char(getAttribute(element,'id'));
name = char(getAttribute(element,'name'));
cancer_name(id) = name;
fprintf('%s\t%s\n',id,name);
n = getLength(subtypes);
subtype_ids = {};
for i=1:n
    nextElement = item(subtypes,i-1);
    if ~isempty(nextElement) && getNodeType(nextElement) == 1
        subid = char(getAttribute(nextElement,'id'));
        subtype_ids = [subtype_ids,subid];
        recurse_type(cancer_tree,cancer_name,nextElement);
    end
    if ~isempty(subtype_ids)
        cancer_tree(id) = subtype_ids;
    end
end
