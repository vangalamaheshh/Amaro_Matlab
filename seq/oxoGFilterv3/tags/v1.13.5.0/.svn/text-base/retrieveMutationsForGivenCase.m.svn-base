function [unameIndices] = retrieveMutationsForGivenCase(uname_case, uname_control, mafTable)
%
% [unameIndices] = retrieveMutationsForGivenCase(uname_case, uname_control, mafTable)
%
% Input:
%   uname_case -- Tumor sample name (usually).  Must match the value of 
%       Tumor_Sample_Barcode inside the maf file.
%
%   uname_control -- Normal sample name (usually).  Must match the value of 
%       Matched_Norm_Sample_Barcode inside the maf file. 
%
%   mafTable -- maf file that should conform to the TCGA MAF Spec.
%               Does require some fields that are not in the MAF spec.
%
%   See also loadMAFTable
%   
% Output:
%   unameIndices -- Indices into the mafTable that correspond to the given
%       case/control combination (columns: Tumor_Sample_Barcode and 
%       Matched_Norm_Sample_Barcode)
%
%
unameIndices = (strcmp(mafTable.Tumor_Sample_Barcode,uname_case) & strcmp(mafTable.Matched_Norm_Sample_Barcode,uname_control));
