function m = convert_varmus_mutlist(infile,outfile,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','hg18');

m = load_struct(infile);

g = load_genes(P.build);

%% patients
m.patient = regexprep(m.Index,'^(.*)$','patient_$1');

% build
m.NCBI_Build = repmat({P.build},slength(m),1);

% coordinates+type
m.chr = convert_chr(m.Chr);
m = make_numeric(m,{'LeftFlank','RightFlank'});
m.start = m.LeftFlank+1;
m.end = m.RightFlank-1;
m.classification = m.muttype;
m.is_indel = strcmp(m.muttype,'INDEL');
m.ref_length = cellfun('length',m.ref_allele);
m.var_length = cellfun('length',m.var_allele);
m.ref_length(strcmp(m.ref_allele,''''''))=0;
m.var_length(strcmp(m.var_allele,''''''))=0;
m.is_ins = m.var_length>m.ref_length;
m.is_del = m.var_length<m.ref_length;
m.classification(m.is_ins) = repmat({'INS'},sum(m.is_ins),1);
m.classification(m.is_del) = repmat({'DEL'},sum(m.is_del),1);
m.is_point = ~m.is_indel & m.ref_length==1 & m.var_length==1;
idx = find(~m.is_point & ~m.is_indel);
if ~isempty(idx)
  fprintf('Now there are apparently DNP''s: please help\n');
  keyboard;
end
m.ref_allele_orig = m.ref_allele;
m.var_allele_orig = m.var_allele;
m.ref_allele = regexprep(m.ref_allele_orig,'''''','-');
m.newbase = regexprep(m.var_allele_orig,'''''','-');
idx = find(m.is_ins);
m.start(idx) = m.LeftFlank(idx)+1;
m.end(idx) = m.RightFlank(idx)-1;
idx = find(m.is_del);
m.gene = m.Gene_name;
idx = find(strcmp('intergenic',m.type));
m.gene(idx) = repmat({'---'},length(idx),1);
idx = grep(',',m.gene,1);

% resolve ambiguous gene names, favoring RefSeq names
for i=1:length(idx)
  names = split(regexprep(regexprep(m.gene{idx(i)},'\"',''),';',','),',');
  is_refseq = ismember(names,g.name);
  if any(is_refseq)
    j = find(is_refseq,1);
    m.gene{idx(i)} = names{j};
  else
    fprintf('gene(s) not found:\n'), disp(names);
  end
end
m.type_orig = m.type;

% convert types
map = {...
    'UTR3','UTR';...
    'UTR3;UTR5','UTR';...
    'UTR5','UTR';...
    'downstream','3''-flank';...
    'frameshift_deletion','frameshift_del';...
    'frameshift_insertion','frameshift_ins';...
    'frameshift_substitution','frameshift';...
    'intergenic','IGR';...
    'intronic','intron';...
    'ncRNA_UTR3','IGR';...
    'ncRNA_UTR5','IGR';...
    'ncRNA_exonic','IGR';...
    'ncRNA_intronic','IGR';...
    'ncRNA_splicing','IGR';...
    'nonframeshift_deletion','inframe_del';...
    'nonframeshift_substitution','inframe';...
    'nonsynonymous_SNV','missense';...
    'splicing','splice_site';...
    'stopgain_SNV','nonsense';...
    'stoploss_SNV','readthrough';...
    'synonymous_SNV','synonymous';...
    'upstream','5''-flank';...
    'upstream;downstream','5''-flank';...
};
aa = []; aa.old = map(:,1); aa.new=map(:,2);
m.type = apply_aliases(m.type_orig,aa);

save_struct(m,outfile);
