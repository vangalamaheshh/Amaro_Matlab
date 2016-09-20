function indel_to_maflite(infile,outfile)

HOMOZYGOUS_CUTOFF = 0.8;

I = load_lines(infile);
I = grepv('REDUCE',I);
I = I(~cellfun('isempty',I));
I = grepv('AUTOFILTER|GERMLINE',I);
I = parse(I,'^(\S+)\t(\S+)\t(\d+)\t(\d+)\t([\-\+])(\S+)\t(.*)$',{'patient','chr','start','stop','insdel','change', 'details'},3:4);

M=[];
M.build = repmat({'36'},slength(I),1);
M.chr=I.chr; M.start=I.start; M.stop=I.stop;
isins = strcmp('+',I.insdel);
M.stop(isins) = M.stop(isins)+1;    % adjust coordinates to the Broad convention
M.start(~isins) = M.start(~isins)+1;
M.ref=I.change; M.ref(isins)=repmat({'-'},sum(isins),1);
M.tum1=I.change; M.tum1(~isins)=repmat({'-'},sum(~isins),1);

cts = parse(I.details,'T_OBS_COUNTS\[C/A/R\]:(\d+)/(\d+)/(\d+)',{'c','a','r'},1:3);
ishom=(cts.c./cts.r)>=HOMOZYGOUS_CUTOFF;
M.tum2=M.ref; M.tum2(ishom)=M.tum1(ishom);

M.tumor_barcode = regexprep(I.patient,'(.*)','$1-Tumor');
M.normal_barcode = regexprep(I.patient,'(.*)','$1-Normal');

save_struct(M,outfile,'no_headers');




return

%% OLD METHOD

P=[];
P.indel_mutation_calls_input_file = infile;
P.indel_mutation_output_maf_file = outfile;
P.input_file_format = 4;
P.output_file_format = 2;
convert_indel_data(P);

