function x = handle_legacy_maflite_format(x,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',true);

% this function only does something if the "official" fieldname is MISSING,
% so actually we don't need to consult P.force_recalculate_maf_simplename_fields

f1 = {'type','classification','start','end','gene','chr','ref_allele','tum_allele1','tum_allele2'};
f2 = {'Variant_Classification','Variant_Type','Start_position','End_position','Hugo_Symbol','Chromosome','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2'};

for i=1:length(f1)
  if isfield(x,f1{i}) && ~isfield(x,f2{i}), x = setfield(x,f2{i},getfield(x,f1{i})); end
end
