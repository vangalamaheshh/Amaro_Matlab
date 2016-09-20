function cn_Dstruct_w_QC = run_CN_QC_analysis(cn_file, QC_output, refgenefile)
% cn_Dstruct_w_QC = run_CN_QC_analysis(cn_file, [QC_output, refgenefile])
% Opens a CN file and reads all data, putting it into a Dstruct object
% Saves QC data to QC_output file
%
% by Michael J.T. O'Kelly, 080416

%% Check inputs

QC_output_default = [cn_file '.QC.test.txt']; %%%
refgenefile_default = '/xchip/tcga/Annotation_Files/UCSC_hg18/hg18_with_miR_20080407.mat';

varlist1 = {'cn_file', 'QC_output', 'refgenefile'};

defaults = {'ERR','QC_output_default','refgenefile_default'};

required = [1,0,0];

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx}])
    end
end

cn_Dstruct = cn_to_Dstruct(cn_file);
load(refgenefile);
cn_Dstruct, cn_Dstruct_w_QC = add_D_snp_scores(cn_Dstruct, cyto);
d = cn_Dstruct_w_QC

% Output QC data

outf = fopen(QC_output, 'w');
fprintf(outf, 'Sample\t%s\t%s\t%s\t%s\t%s\n', d.supdesc(1,:), d.supdesc(2,:), d.supdesc(3,:), d.supdesc(4,:), d.supdesc(5,:))
supdat_size = size(d.supdat)
for row_index = 1:supdat_size(2)
	sample_name = cell2mat(d.sdesc(row_index))
	row_data = d.supdat(:,row_index)
	fprintf(outf, '%s\t%f\t%f\t%f\t%f\t%f\n', sample_name, row_data(1), row_data(2), row_data(3), row_data(4), row_data(5))
end

fclose(outf);


