function Y = convert_mutation_file_1(X,P)
if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'tumor_type_prefix','OV');
P=impose_default_value(P,'dataset','CAP');

Y = [];
Y.chr = regexprep(X.Chromosome,'(.*)','chr$1');
Y.start = X.Start_position;
Y.end = X.End_position;
Y = make_numeric(Y,{'start','end'});
Y.patient = regexprep(X.Tumor_Sample_Barcode,'^TCGA-..-(....)-.*$',[P.tumor_type_prefix '-$1']);
a = {'Missense_Mutation','Nonsense_Mutation',...
   'frame_shift_del','frame_shift_ins','in_frame_del','in_frame_ins',...
   'nonstop','silent','splice_site','splice_site_del','splice_site_ins','Splice_Site_SNP',...
   'Splice_Site_Indel','Frame_Shift_Del','Frame_Shift_Ins','Silent',...
};
b = {'Missense','Nonsense',...
   'Frameshift_Del','Frameshift_Ins','Inframe_Del','Inframe_Ins',...
   'Nonsense','Synonymous','Splice_site','Frameshift_Del','Frameshift_Ins','Splice_site',...
   'Splice_site','Frameshift_Del','Frameshift_Ins','Synonymous',...
};

Y.type = map_across(X.Variant_Classification,a,b);
Y.class = upper(X.Reference_Allele);
idx = union(grep('In|Del',Y.type,1),find(Y.end-Y.start>0));
Y.class(idx) = repmat({'Indel'},length(idx),1);
% annotate "class" (base categories)
for i=1:slength(Y)
  if strcmp(Y.class{i},'C')
    base = upper(genome_region(Y.chr{i},Y.start(i)+1));
    if strcmp(base,'G'), Y.class{i} = 'C (CpG)';
    else Y.class{i} = 'C (other)'; end
  elseif strcmp(Y.class{i},'G')
    base = upper(genome_region(Y.chr{i},Y.start(i)-1));
    if strcmp(base,'C'), Y.class{i} = 'G (CpG)';
    else Y.class{i} = 'G (other)'; end
  end
end
Y.site = cell(slength(X),1);
for i=1:slength(X), Y.site{i} = [X.Hugo_Symbol{i} '|chr' X.Chromosome{i} ':' X.Start_position{i} '-' X.End_position{i}]; end
Y.gene = X.Hugo_Symbol;
Y.ref_allele = X.Reference_Allele;
Y.tum_allele1 = X.Tumor_Seq_Allele1;
Y.tum_allele2 = X.Tumor_Seq_Allele2;
Y.transcript = repmat({'-'},slength(X),1);
Y.proteinchange = repmat({'-'},slength(X),1);
Y.dataset = repmat({P.dataset},slength(X),1);
Y.filtered = repmat({'OK'},slength(X),1);
