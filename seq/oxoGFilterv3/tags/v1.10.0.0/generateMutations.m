function [R, R_A, R_nA] = generateMutations(mafTable, unameIndices)
%
% [R, R_A, R_nA] = generateMutations(mafTable, unameIndices)
% R -- combined list of R_A and R_nA with order preserved from the inpuit
%   mafTable
% R_A -- indices of artifact mode mutations.  Note that this will include real
%   mutations and artifacts.  Simply a collection of C>A/G>T mutations.
% R_nA -- indices of mutations that are not in the artifact mode
% Generate mutations as a single struct.
%
% Each field in the struct is a column vector of values with length being 
%   equal to the number of mutations.  Each field has the same length and 
%   each row corresponds to the same mutation across fields.
%
%   The struct contains the following fields:
%
%   foxog -- proportion of alt reads in the oxoG orientation
%   alt_read_count -- total number of alternate reads
%   ref_read_count -- total number of reference reads
%   f1r2 -- number of alt reads with f1r2 orientation
%   f2r1 -- number of alt reads with f2r1 orientation
%   isArtifactMode -- true/false (1/0) whether this mutation qualifies as
%       artifact mode (i.e. C>A or G>T)
%
% f1r2 + f2r1 will be equal to alt_read_count
% round(foxog * alt_read_count) will give an integer number of the count of
%   oxoG orientation reads in the alt reads.  
%
% Command to count artifact mode mutations for a given barcode (e.g.
%   HSCX1127)
% egrep HSCX1127 PR_SIGMA_Cervical_Capture.snp.oxoG.tcga.all.maf.annotated | egrep "SNP\WC\WA|SNP\WG\WT" | wc
%
% Note:  This function creates two copies of the mutations.
%  That may be an issue in some cases due to RAM consumption.
%

refAlleles = mafTable.Reference_Allele(unameIndices);
altAlleles = mafTable.Tumor_Seq_Allele1(unameIndices);
unameFoxoG = mafTable.i_t_Foxog(unameIndices);
unameAltF1R2 = mafTable.i_t_ALT_F1R2(unameIndices);
unameAltF2R1 = mafTable.i_t_ALT_F2R1(unameIndices);
unameRefF1R2 = mafTable.i_t_REF_F1R2(unameIndices);
unameRefF2R1 = mafTable.i_t_REF_F1R2(unameIndices);


R.foxog = unameFoxoG;
R.alt_read_count = (unameAltF1R2+unameAltF2R1);
R.ref_read_count = (unameRefF1R2+unameRefF2R1);
R.f1r2 = unameAltF1R2;
R.f2r1 = unameAltF2R1;
R.isArtifactMode = isArtifactSignature(refAlleles, altAlleles);

% As part of maf file convention, non-artifact were being populated with -1
%   in old versions, so recalculate any foxog of -1.
% 
% To build null-model FoxoG, choose an arbitrary read config (e.g.
%   F2R1) and reconstruct a FoxoG.  Any biases here are being ignored.
%
% Convention 
%   C or A to anything: F2R1
%   G or T to anything: F1R2

R.foxog = unameFoxoG;
f2r1_i = (strcmp(refAlleles, 'C') | strcmp(refAlleles, 'A')) & (R.foxog == -1);
f1r2_i = (strcmp(refAlleles, 'G') | strcmp(refAlleles, 'T')) & (R.foxog == -1);
R.foxog(f2r1_i) = R.f2r1(f2r1_i) ./ R.alt_read_count(f2r1_i);
R.foxog(f1r2_i) = R.f1r2(f1r2_i) ./ R.alt_read_count(f1r2_i);

R.alt_read_count = (unameAltF1R2+unameAltF2R1);
R.ref_read_count = (unameRefF1R2+unameRefF2R1);
R.f1r2 = unameAltF1R2;
R.f2r1 = unameAltF2R1;

% Create R_A and R_nA
R_A_i = isArtifactSignature(refAlleles, altAlleles);
R_nA_i = ~isArtifactSignature(refAlleles, altAlleles);

R_A.foxog = R.foxog(R_A_i);
R_A.alt_read_count = (unameAltF1R2(R_A_i)+unameAltF2R1(R_A_i));
R_A.ref_read_count = (unameRefF1R2(R_A_i)+unameRefF2R1(R_A_i));
R_A.f1r2 = unameAltF1R2(R_A_i);
R_A.f2r1 = unameAltF2R1(R_A_i);
R_A.isArtifactMode = R.isArtifactMode(R_A_i);
if any(~R_A.isArtifactMode)
   error('Bug detected.  Artifact mode mutations should all be in the artifact mode') 
end


R_nA.foxog = R.foxog(R_nA_i);
R_nA.alt_read_count = (unameAltF1R2(R_nA_i)+unameAltF2R1(R_nA_i));
R_nA.ref_read_count = (unameRefF1R2(R_nA_i)+unameRefF2R1(R_nA_i));
R_nA.f1r2 = unameAltF1R2(R_nA_i);
R_nA.f2r1 = unameAltF2R1(R_nA_i);
R_nA.isArtifactMode = R.isArtifactMode(R_nA_i);
if any(R_nA.isArtifactMode)
   error('Bug detected.  Non-Artifact mode mutations should not be in the artifact mode') 
end

