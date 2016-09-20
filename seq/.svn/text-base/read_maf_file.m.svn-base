function maf=read_maf_file(fname,output_format)
%
% read_maf_file(fname,output_format)
%
% Gaddy format:
% if output_format is omitted, returns maf as struct with fields
%    sdesc, gacc, gdesc, gsymb, gdat
%
% Mike format:
% if output_format is "struct", returns maf as struct with fields
%    matching the columns of the maf file
%

if ~exist('output_format', 'var')
  output_format = '';
end

if ~strcmp(output_format, 'struct')

  tbl=read_table(fname,[],char(9),1,'Whitespace',' \b\r');

  flds={'Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_position','End_position', ...
      'Strand','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele1',...
      'Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode',...
      'Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2',...
      'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status',...
      'Mutation_Status'};

  [Mt,m1,m2]=match_string_sets_hash(lower(tbl.headers{1}),lower(flds));
  if ~isunique(m2) || length(m2)~=length(flds)
    analyze_string_sets(lower(tbl.headers{1}),lower(flds),[1 3]);
  %  keyboard
  end

  % keyboard

  maf.sdesc=flds;
  hugo_idx=find(m2==1);
  if isempty(hugo_idx)
    keyboard
    error('Cant find Hugo_Symbol');
  end
  hugo_idx=m1(hugo_idx);
  maf.gacc=tbl.dat{hugo_idx(1)};
  maf.gdesc=maf.gacc;
  maf.gsymb=maf.gacc;
  maf.dat=cell(length(maf.gacc),length(flds));

  for i=1:length(flds)
    idx=find(m2==i);
    if isempty(idx)
      maf.dat(:,i)=cellstr(repmat('NULL',length(maf.gacc),1));
      verbose(['filling column ' flds{i} ' with NULL'],30);
    else
      maf.dat(:,i)=tbl.dat{m1(idx)};
    end
  end

else    % output_format == 'struct'

  field = {       % mapping from structure fields to column names

    'patient', 'tumor_sample' ;...
    'gene', 'hugo_symbol' ;...
    'chr', 'chrom' ;...
    'start', 'start_position' ;...
    'end', 'end_position' ;...
    'type', 'variant_classification' ;...
    'status', 'mutation_status' ;...
    'verification', 'verification_status' ;...
    'validation', 'validation_status' ;...
    'center', 'center' ;...

    'ref', 'reference_allele' ;...
    'ts1', 'tumor_seq_allele1' ;...
    'ts2', 'tumor_seq_allele2' ;...
    'mns1', 'Match_Norm_Seq_Allele1' ;...
    'mns2', 'Match_Norm_Seq_Allele2' ;...
    'tv1', 'Tumor_Validation_Allele1' ;...
    'tv2', 'Tumor_Validation_Allele2' ;...
    'mnv1', 'Match_Norm_Validation_Allele1' ;...
    'mnv2', 'Match_Norm_Validation_Allele2' ;...
  };

  table = read_table(fname,repmat('%s',1,28),char(9),1,'whitespace','\b\r');
  maf = [];
  for f=1:size(field,1)
    col = find(strncmpi(table.headers{1},field{f,2},length(field{f,2})));
    if isempty(col), error('No %s column in %s\n', field{f,2}, fname); end
    maf = setfield(maf, field{f,1}, table.dat{col});
  end

end

