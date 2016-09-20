function x = add_and_convert_simple_fieldnames(x,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',false);

try
  x = add_simple_fieldnames(x,P);
catch me
  disp(me.message);
  fprintf('Warning: add_simple_fieldnames failed.  Maybe fieldnames are already simple.\n');
end

if isfield(x,'patient_name') && isfield(x,'patient'),
  % if patient looks like it's really patient_idx, then replace it with patient_name
  tmp = str2double(x.patient);
  if mean(isnan(tmp))<0.05
    x.patient = x.patient_name;
  end
end

if isfield(x,'gene_name') && isfield(x,'gene'),
  tmp = str2double(x.gene);
  if mean(isnan(tmp))<0.05
    x.gene = x.gene_name;
  end
end

if ~isfield(x,'newbase'), x.newbase = find_newbase(x); end

mandatory_flds = {'patient','chr','start','end','type','gene','classification','ref_allele','newbase'};
require_fields(x,mandatory_flds);

x.chr = convert_chr(x.chr,P);
x = make_numeric(x,{'start','end'});

if isfield(x,'firehose_patient_id'), x.patient = x.firehose_patient_id; end   % MAKE SURE



