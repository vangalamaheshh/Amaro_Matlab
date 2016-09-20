function indel_to_tcgamaf(infile,outfile,seqsource)

if ~exist('seqsource','var'), seqsource = '---'; end

HOMOZYGOUS_CUTOFF = 0.8;

I = load_lines(infile);
I = grepv('REDUCE',I);
I = I(~cellfun('isempty',I));
I = grepv('AUTOFILTER|GERMLINE',I);
I = parse(I,['^(\S+)\t(\S+)\t(\d+)\t(\d+)\t([\-\+])(\S+)\t(' repmat('\S+\t',1,11) '\S+)\t(\S+)\t(\S+)\t?(\S*)$'],...
  {'patient','chr','start','stop','insdel','change','details','somatic','type','gene'},3:4);

M=[];
M.Hugo_Symbol = fillblanks(I.gene,'---');
z = repmat({'---'},slength(M),1);
M.Entrez_Gene_Id = z;
M.Center = repmat({'broadinstitute.org'},slength(M),1);
M.NCBI_Build = repmat({'36'},slength(I),1);
M.Chromosome = regexprep(I.chr,'^chr','');
M.Start_position = I.start; 
M.End_position = I.stop;
M.Strand = repmat({'+'},slength(M),1);

isinframe = mod(cellfun('length',I.change),3)==0;
isins = strcmp('+',I.insdel);
M.End_position(isins) = M.End_position(isins)+1;    % adjust coordinates to the Broad convention
M.Start_position(~isins) = M.Start_position(~isins)+1;

M.Variant_Classification = stringsplice([nansub({'Frame_Shift';'In_Frame'},isinframe+1),...
  nansub({'Del';'Ins'},isins+1)],1,'_');
idx = grep('GENOMIC',I.type,1); M.Variant_Classification(idx) = repmat({'IGR'},length(idx),1);
idx = grep('INTRON',I.type,1); M.Variant_Classification(idx) = repmat({'Intron'},length(idx),1);
idx = grep('UTR|UNKNOWN',I.type,1); M.Variant_Classification(idx) = repmat({'UTR'},length(idx),1);
M.Variant_Type = nansub({'DEL';'INS'},isins+1);

M.Reference_Allele = I.change; M.Reference_Allele(isins)=repmat({'-'},sum(isins),1);
M.Tumor_Seq_Allele1 = I.change; M.Tumor_Seq_Allele1(~isins)=repmat({'-'},sum(~isins),1);
cts = parse(I.details,'T_OBS_COUNTS\[C/A/R\]:(\d+)/(\d+)/(\d+)',{'c','a','r'},1:3);
ishom = (cts.c./cts.r)>=HOMOZYGOUS_CUTOFF;
M.Tumor_Seq_Allele2 = M.Reference_Allele; M.Tumor_Seq_Allele2(ishom) = M.Tumor_Seq_Allele1(ishom);

M.dbSNP_RS = z;
M.dbSNP_Val_Status = z;

M.Tumor_Sample_Barcode = regexprep(I.patient,'(.*)','$1-Tumor');
M.Matched_Norm_Sample_Barcode = regexprep(I.patient,'(.*)','$1-Normal');

M.Match_Norm_Seq_Allele1 = z;
M.Match_Norm_Seq_Allele2 = z;
M.Tumor_Validation_Allele1 = z;
M.Tumor_Validation_Allele2 = z;
M.Match_Norm_Validation_Allele1 = z;
M.Match_Norm_Validation_Allele2 = z;
M.Verification_Status = z;
M.Validation_Status = z;
M.Mutation_Status = repmat({'Somatic'},slength(M),1);
M.Sequencing_Phase = z;
M.Sequence_Source = repmat({seqsource},slength(M),1);
M.Validation_Method = z;
M.Score = z;
M.BAM_file = z;
M.Sequencer = z;

save_struct(M,outfile);

