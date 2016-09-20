function INDEL = load_indel_data(P)
% reads indel call file
%
% returns struct INDEL with the following fields:
%  patient
%  chr,start,end
%  ref_allele
%  tum_allele1
%  tum_allel2
%  type
%  gene
%  dataset
%
% Mike Lawrence 2009-07-02

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'capture_patient_list_input_file','*required*');
P=impose_default_value(P,'indel_mutation_list_input_file','*required*');
P=impose_default_value(P,'indel_homozygous_cutoff_ratio',0.7);

CAP_patients = load_struct(P.capture_patient_list_input_file);

INDEL = load_struct_specify_numeric_cols(P.indel_mutation_list_input_file,[],0);
INDEL = rename_fields(INDEL,colx(1:8),{'file','chr','start','end','details','somatic_status','coding_status','gene'});
INDEL = reorder_struct(INDEL,strcmp('CODING',INDEL.coding_status));
INDEL = reorder_struct(INDEL,strcmp('SOMATIC',INDEL.somatic_status));
tmp = parse(INDEL.details,'^(-|\+)([^:]*):(\d*)/(\d*)$',{'sign','seq','affected_reads','total_reads'});
tmp = make_numeric(tmp,{'affected_reads','total_reads'});
INDEL.is_insertion = strcmp(tmp.sign,'+');
INDEL.is_inframe = (mod(cellfun('length',tmp.seq),3)==0);
INDEL.is_homozygous = ((tmp.affected_reads ./ tmp.total_reads) > P.indel_homozygous_cutoff_ratio);
indel_type = {'Frameshift_Del','Frameshift_Ins','Inframe_Del','Inframe_Ins'};
INDEL.type = indel_type(1+(INDEL.is_insertion + 2*INDEL.is_inframe))';
short = regexprep(INDEL.file,'^.*/OV-(\d\d\d\d)/.*$','$1');
idx = listmap(short,CAP_patients.short);
if any(isnan(idx))
  fprintf('The following patients were not found in %s\n',P.CAP_patient_list);
  disp(short(isnan(idx)));
  error('Please add them and try re-running the analysis\n');
end
INDEL.patient = CAP_patients.name(idx);
INDEL.chr = convert_chr(INDEL.chr);
INDEL.dataset = repmat({'CAP'},slength(INDEL),1);
INDEL.ref_allele = repmat({'-'},slength(INDEL),1);
INDEL.tum_allele1 = repmat({'-'},slength(INDEL),1);
INDEL.tum_allele2 = repmat({'-'},slength(INDEL),1);
idx = find(INDEL.is_insertion);
INDEL.tum_allele1(idx) = tmp.seq(idx);
idx = intersect(idx,find(INDEL.is_homozygous));
INDEL.tum_allele2(idx) = tmp.seq(idx);
idx = find(~INDEL.is_insertion);
INDEL.ref_allele(idx) = tmp.seq(idx);
idx = intersect(idx,find(~INDEL.is_homozygous));
INDEL.tum_allele1(idx) = tmp.seq(idx);
