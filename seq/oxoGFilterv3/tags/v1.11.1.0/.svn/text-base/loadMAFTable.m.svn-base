function [mafTable] = loadMAFTable(mafFilename)
%
% [mafTable] = loafMAFTable(mafFilename)
%
% Loads a maf file as a table.  Structure is the same as load_table.
%
% Note:  This loads the entire table.  However, only certain columns are required.
%   Those are:
%
%   Chromosome
%   Start_position
%   End_position
%   Variant_Type (only for lego plots)
%   Reference_Allele
%   Tumor_Seq_Allele1
%   Tumor_Sample_Barcode
%   Matched_Norm_Sample_Barcode
%   validation_alt_allele
%   i_t_ALT_F1R2
%   i_t_ALT_F2R1
%   i_t_REF_F1R2
%   i_t_REF_F2R1
%   i_t_Foxog
%   ref_context (only for lego plots)
%
% These correspond to columns: [5:7 10 11 12 16 17 64 83 90:95]
%   
% Note:  No field can have a value larger than 16k in a single entry
%
% See also load_table
%
if ~exist(mafFilename, 'file')
   error(['Could not find maf file: ' mafFilename])
end

% Need to determine the number of header lines as this can change and 2 is
%   currently hardcoded.
mafTable = load_table(mafFilename, char(9));

% Sanity check:
requiredHeaders = [{'Chromosome'}; {'Start_position'}; {'End_position'}; {'Reference_Allele'}; {'Tumor_Seq_Allele1'}; {'Tumor_Sample_Barcode'}; {'Matched_Norm_Sample_Barcode'}; {'ref_context'}; {'validation_alt_allele'}; {'i_t_ALT_F1R2'}; {'i_t_ALT_F2R1'}; {'i_t_REF_F1R2'}; {'i_t_REF_F2R1'}; {'i_t_Foxog'}];
for i = 1:length(requiredHeaders)
   requiredHeader = requiredHeaders{i};
   check = sum( strcmp(mafTable.header, requiredHeader) );
   if check == 0
      error(['Unable to find required header: '  requiredHeader '      Found: ' mafTable.header '      Required Headers: ' requiredHeaders])
   end
end
