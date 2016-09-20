function convert_indel_data(P)
% if P.input_file_format == 1
%
%  reads indel call file with the following columns (no header row):
%    (assumes all are CODING and SOMATIC)
%     patient
%     chr
%     start
%     end
%     details (e.g. "+CTT:562/592")
%     gene
%     [verification_status]
%     [validation_status]
%
% if P.input_file_format == 2
%
%  reads indel call file with the following columns (no header row):
%    (assumes all are CODING and SOMATIC)
%     patient
%     chr
%     start
%     end
%     details (e.g. "+CTT:562/592")
%     12 columns of stats
%     [somatic status] --> filter to include only "SOMATIC"
%     [coding status]  --> filter to include only "CODING"
%     gene
%
% if P.input_file_format == 3
%
%  reads indel call file with the following columns (no header row):
%    (assumes all are CODING and SOMATIC)
%     patient
%     chr
%     start
%     end
%     details (e.g. "+CTT:562/592")
%     12 columns of stats
%     [somatic status] --> filter to include only "SOMATIC"
%     [coding status]  --> filter to include only "CODING"
%     gene
%     autofilter
%
% writes indel maf file with the following columns (with header row):
%   patient
%   chr
%   start
%   end
%   ref_allele
%   tum_allele1
%   tum_allel2
%   type {'Frameshift_Del','Frameshift_Ins','Inframe_Del','Inframe_Ins'}
%   gene
%   [verification_status]
%   [validation_status]
%
% Mike Lawrence 2009-11-09

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'indel_mutation_calls_input_file','*required*');
P=impose_default_value(P,'indel_mutation_output_maf_file','*required*');
P=impose_default_value(P,'indel_homozygous_cutoff_ratio',0.7);
P=impose_default_value(P,'input_file_format',1);
P=impose_default_value(P,'output_file_format',1);
P=impose_default_value(P,'dataset','INDEL');
P=impose_default_value(P,'build','hg18');

if P.input_file_format==1
  I = load_struct_noheader(P.indel_mutation_calls_input_file);
  I = rename_fields(I,colx(1:6),{'patient','chr','start','end','details','gene'});
  if isfield(I,'col7'), I = rename_field(I,'col7','verification_status');
  else I.verification_status = repmat({'UNKNOWN'},slength(I),1); end
  if isfield(I,'col8'), I = rename_field(I,'col8','validation_status');
  else I.validation_status = repmat({'UNKNOWN'},slength(I),1); end
elseif P.input_file_format==2
  I = load_struct_noheader(P.indel_mutation_calls_input_file);
  I = rename_fields(I,colx(1:5),{'patient','chr','start','end','details'});
  I = rename_fields(I,colx(18:20),{'somatic_status','coding_status','gene'});
  I = reorder_struct(I,strcmp('SOMATIC',I.somatic_status));
  I = reorder_struct(I,strcmp('CODING',I.coding_status));
  I.verification_status = repmat({'UNKNOWN'},slength(I),1);
  I.validation_status = repmat({'UNKNOWN'},slength(I),1);
elseif P.input_file_format==3
  I = load_struct_noheader(P.indel_mutation_calls_input_file,24);
  I = rename_fields(I,colx(1:5),{'patient','chr','start','end','details'});
  I = rename_fields(I,colx(18:24),{'somatic_status','coding_status','gene','filter',...
     'filter2','filter3','filter4'});
  I = reorder_struct(I,strcmp('SOMATIC',I.somatic_status));
  I = reorder_struct(I,strcmp('CODING',I.coding_status));
  idx = grep('^AUTOFILTER|FAILED_REVIEW',I.filter,1);
  I = reorder_struct(I,setdiff(1:slength(I),idx));
  I.verification_status = repmat({'UNKNOWN'},slength(I),1);
  I.validation_status = repmat({'UNKNOWN'},slength(I),1);
  I.filters = stringsplice([I.filter2,I.filter3,I.filter4],1);
  idx = grep('REVIEW',I.filter,1);
  I.verification_status(idx) = I.filter(idx);
  I.filters = stringsplice([I.filter2 I.filter3 I.filter4]);
  idx = grep('SQNM',I.filters,1);
  I.validation_status(idx) = I.filters(idx);
elseif P.input_file_format==4
  I = load_struct_noheader(P.indel_mutation_calls_input_file,24);
  I = rename_fields(I,colx(1:5),{'patient','chr','start','end','details'});
  I = rename_fields(I,colx(18:24),{'somatic_status','coding_status','gene','filter',...
     'filter2','filter3','filter4'});
  idx = grep('^AUTOFILTER|FAILED_REVIEW',I.filter,1);
  I = reorder_struct(I,setdiff(1:slength(I),idx));
  I.verification_status = repmat({'UNKNOWN'},slength(I),1);
  I.validation_status = repmat({'UNKNOWN'},slength(I),1);
  I.filters = stringsplice([I.filter2,I.filter3,I.filter4],1);
  idx = grep('REVIEW',I.filter,1);
  I.verification_status(idx) = I.filter(idx);
  I.filters = stringsplice([I.filter2 I.filter3 I.filter4]);
  idx = grep('SQNM',I.filters,1);
  I.validation_status(idx) = I.filters(idx);
else
  error('unknown P.input_file_format');
end

for i=1:slength(I), if ~contains(I.details{i},':'), I.details{i} = [I.details{i} ':1/2']; end, end
tmp = parse(I.details,'^(-|\+)([^:]*):(\d*)/(\d*)$',{'sign','seq','affected_reads','total_reads'});
tmp = make_numeric(tmp,{'affected_reads','total_reads'});
I.is_insertion = strcmp(tmp.sign,'+');
I.is_inframe = (mod(cellfun('length',tmp.seq),3)==0);
I.is_homozygous = ((tmp.affected_reads ./ tmp.total_reads) > P.indel_homozygous_cutoff_ratio);
indel_type = {'Frameshift_Del','Frameshift_Ins','Inframe_Del','Inframe_Ins'};
I.type = indel_type(1+(I.is_insertion + 2*I.is_inframe))';
I.chr = convert_chr(I.chr);
I = reorder_struct(I,~isnan(I.chr));
I.ref_allele = repmat({'-'},slength(I),1);
I.tum_allele1 = repmat({'-'},slength(I),1);
I.tum_allele2 = repmat({'-'},slength(I),1);
idx = find(I.is_insertion);
I.tum_allele1(idx) = tmp.seq(idx);
idx = intersect(idx,find(I.is_homozygous));
I.tum_allele2(idx) = tmp.seq(idx);
idx = find(~I.is_insertion);
I.ref_allele(idx) = tmp.seq(idx);
idx = intersect(idx,find(~I.is_homozygous));
I.tum_allele1(idx) = tmp.seq(idx);
from = {'NOT_REVIEWED','FAILED','GERMLINE','VALID','POSSIBLY_VALID'};
to = {'NOT_REVIEWED','FAILED_REVIEW','FAILED_REVIEW(seen in normal)','PASSED_REVIEW','PASSED_REVIEW(weak)'};
for i=1:length(from)
  idx = find(strcmp(I.verification_status,from{i}));
  I.verification_status(idx) = repmat(to(i),length(idx),1);
end

I.build = repmat({'36'},slength(I),1);
I.chr = chrlabel(I.chr);
I.classification = repmat({'Indel'},slength(I),1);

if P.output_file_format == 1;
  I.dataset = repmat({P.dataset},slength(I),1);
  flds = {'patient','build','chr','start','end','ref_allele','tum_allele1','tum_allele2','type',...
  'classification','gene','dataset','verification_status','validation_status'};
  I = keep_fields(I,flds);
  I = orderfields(I,flds);
elseif P.output_file_format == 2;
  [I.tumor_barcode I.normal_barcode] = get_TCGA_barcodes(I.patient);
  I.validation_method = repmat({''},slength(I),1);
  idx = grep('SQNM',I.validation_status,1);I.validation_method(idx) = repmat({'Sequenom'},length(idx),1);
  idx1 = grep('SQNM_SOMATIC',I.validation_status,1);I.validation_status(idx1) = repmat({'Validated'},length(idx1),1);
  idx2 = grep('SQNM_GERMLINE',I.validation_status,1);I.validation_status(idx2) = repmat({'Germline'},length(idx2),1);
  idx3 = setdiff(1:slength(I),[idx1;idx2]); I.validation_status(idx3) = repmat({'Unknown'},length(idx3),1);
  I.verification_status = repmat({''},slength(I),1);
  flds = {'build','chr','start','end','ref_allele','tum_allele1','tum_allele2',...
          'tumor_barcode','normal_barcode','gene','type','classification','verification_status','validation_status'};
  I = keep_fields(I,flds);
  I = orderfields(I,flds);
  % need to get transcript from Refseq (so Entrez gene ID can be loaded in Alex's script)
  R = load_refseq(P.build);
  I.transcript = repmat({''},slength(I),1);
  for i=1:slength(I)
    idx = find(strcmpi(R.gene,I.gene{i}),1);
    if ~isempty(idx), I.transcript{i} = R.transcript{idx}; end
  end
else
  error('Unknown P.output_file_format');
end

save_struct(I,P.indel_mutation_output_maf_file);
