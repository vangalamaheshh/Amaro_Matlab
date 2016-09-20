function varargout = check_for_absolute_segs(input)
%% Checks whether the arrays have been run through absolute.
% Two possible inputs: 
% 1. Cell array of SNP array names. 
% 2. SIF file path with two columns: Individual ID and SNP array name
% Checks in the following directory: ~scarter/bioinf/birdseed_PL/proc_seg/ABSOLUTE_db/d1_CEPH_MN_IMP_Mseg_BG_pBSv0.6_H1_t/CALLED/SNP_6.0]
%
% If you provided just snp array names returns the absolute path to the ABSOLUTE seg file and the R object. 
% If you provided a path to a SIF file, returns the individual names as the third argument. 
% Petar Stojanov - 11/30/2011

individuals = false;
if ischar(input)
  try
    sif = load_struct_noheader(input);
    array_list = sif.col2;
    individual_ids = sif.col1;
    individuals = true;
  catch
    error('Problem with path to SIF');
  end
elseif iscell(input)
  array_list = input;
else 
  error('Must provide either path to SIF or cell array of SNP array names');
end


plates = cellfun(@splitstr, array_list, repmat({'_'}, length(array_list), 1), 'UniformOutput', false);
plates = cellfun(@(x) x{1}, plates, 'UniformOutput', false);
absolute_paths = cellfun(@(x, y) ['~scarter/bioinf/birdseed_PL/proc_seg/ABSOLUTE_db/d1_CEPH_MN_IMP_Mseg_BG_pBSv0.6_H1_t/CALLED/SNP_6.0/' ...
                    x '/' y '__ABSOLUTE_CALLED_segtab.txt'], plates, array_list,'UniformOutput', false);
absolute_paths_R = cellfun(@(x, y) ['~scarter/bioinf/birdseed_PL/proc_seg/ABSOLUTE_db/d1_CEPH_MN_IMP_Mseg_BG_pBSv0.6_H1_t/CALLED/SNP_6.0/' ...
                    x '/' y '__ABSOLUTE_CALLED.Rdata'], plates, array_list,'UniformOutput', false);


exists = cellfun(@fopen, absolute_paths, 'UniformOutput', false);
exists = cellfun(@double, exists);
idx = find(exists ~= -1);
varargout{1} = absolute_paths(idx);
varargout{2} = absolute_paths_R(idx);
if individuals
  varargout{3} = individual_ids(idx);
end