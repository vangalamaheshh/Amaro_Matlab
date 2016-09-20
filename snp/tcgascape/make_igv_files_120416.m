% make igv files for 120416 run

dbstop if error
set_verbose_level(30);

% method for determining disease type
%   true  = use a sample info on single pan-cancer seg file
%   false = use separate segfiles for each disease, mapped by
%           tcga_segfgiles()
sample_info_method = false;

output_dir = '/xchip/gistic/tcgascape/tcgascape_120416/';

% cache file for cleaned D
cached_clean_master = [output_dir,'Dmaster.mat'];

% create output directory
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% load cancer type hierarchy and name mappings
[cancer_namemap cancer_treegen] = tcga_cancer_types();

D = load_D(cached_clean_master);

%% save D structs and segfiles for IGV
segfile_dir = [output_dir 'igvfiles/'];

if ~exist(segfile_dir,'dir');
    mkdir(segfile_dir);
    unix(['chmod 775 ' segfile_dir]);
end

write_igv_segs(D,'all_cancers',cancer_treegen,segfile_dir,cancer_namemap,40);
