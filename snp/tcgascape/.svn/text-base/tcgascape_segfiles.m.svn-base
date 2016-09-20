function input_paths = tcgascape_segfiles()
%TCGASCAPE_SEGFILES map TCGA-define cancer types to segfile input paths
%
%   INPUT_PATHS = TCGASCAPE_SEGFILES()
%
% Returns a containers.Map from TCGA-defined tissue types to the segmented
% copy number data input files for each type. This function should be 
% modified when the input data files change (i.e. every run).

% 21-Mar-2012 run - use programmatically named seg files from Nam

% segment directory, copied from /xchip/cga1/nampho/seg_v2.1 3/21/2012 6:00 PM
seg_dir = '/xchip/gistic/tcgascape/tcgascape_120321/input_segs/';
seg_files = dir(seg_dir);

fnames = regexp({seg_files.name},'.+_hg19_T_qc.merged.seg','match');
fnames = fnames(~cellfun(@isempty,fnames));
fnames = [fnames{:}];
types = regexprep(fnames,'(.+)_hg19_T_qc.merged.seg','$1');
input_paths = containers.Map(types,strcat(seg_dir,fnames));
