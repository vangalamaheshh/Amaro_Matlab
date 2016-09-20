function x = move_to_simple_fieldnames(x,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',true);

flds1 = {'Variant_Classification','Variant_Type','Start_position','End_position','Hugo_Symbol','Chromosome','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2'};
flds2 = {'type','classification','start','end','gene','chr','ref_allele','tum_allele1','tum_allele2'};

for i=1:length(flds1)
  if isfield(x,flds1{i})
    if isfield(x,flds2{i})
      tmp = getfield(x,flds1{i});
      empty1 = false;
      if iscell(tmp) && all(strcmp('',tmp)) empty1 = true; end;
      if isnumeric(tmp) && all(isnan(tmp)), empty1 = true; end;
      empty2 = false;
      if iscell(tmp) && all(strcmp('',tmp)) empty2 = true; end;
      if isnumeric(tmp) && all(isnan(tmp)), empty2 = true; end;
      if P.force_recalculate_maf_simplename_fields, empty2 = true; end
      if empty1
        x = rmfield(x,flds1{i});
      elseif empty2
        x = rmfield(x,flds2{i});
        x = rename_field(x,flds1{i},flds2{i});
      else
        fprintf('Warning has both %s and %s; overwriting %s with %s\n',flds1{i},flds2{i},flds2{i},flds1{i});
        x = rmfield(x,flds2{i});
        x = rename_field(x,flds1{i},flds2{i});
      end
    else
      x = rename_field(x,flds1{i},flds2{i});
    end
  end
end

