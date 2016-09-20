function p = dRanger_design_primers_from_validation_targets(P)
% Mike Lawrence 2009-10-26

if ~exist('P','var'), P = []; end

P = impose_default_value(P,'outstem','*required*');
P = impose_default_value(P,'include_A_and_B_controls',false);
P = impose_default_value(P,'mask_central_region',150);

if P.include_A_and_B_controls, error('P.include_A_and_B_controls not currently supported.'); end

rfile = [P.outstem '_rearrs.txt'];
tfile = [P.outstem '_T.fasta'];
demand_file({rfile,tfile});

T = load_fasta(tfile);
R = load_struct(rfile);
require_fields(R,{'seqname'});

seqlen = cellfun('length',T.seq);
if any(seqlen<P.mask_central_region), error('Some sequences are shorter than P.mask_central_region'); end

primer3_params = {...
   'PRIMER_TASK=pick_pcr_primers';   % 'pick_sequencing_primers' doesn't work at all!
   'PRIMER_FIRST_BASE_INDEX=1';
   'PRIMER_PICK_LEFT_PRIMER=1';
   'PRIMER_PICK_RIGHT_PRIMER=1';
   'PRIMER_NUM_RETURN=1';
   'PRIMER_MIN_SIZE=18';
   'PRIMER_OPT_SIZE=21';
   'PRIMER_MAX_SIZE=27';
   'PRIMER_PRODUCT_SIZE_RANGE=250-400';
   'PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.05';
   'PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.05';
   'PRIMER_PRODUCT_OPT_SIZE=325';
   'PRIMER_MAX_NS_ACCEPTED=0';
   'PRIMER_MIN_GC=40';
   'PRIMER_MAX_GC=60';
   'PRIMER_MAX_POLY_X=10';  
};

targleft = round((seqlen - P.mask_central_region)/2);
seq_input = cell(slength(T),1);
for i=1:slength(T)
  name = T.header{i}; name(name==13) = [];
  seq = T.seq{i}; seq(seq==13) = [];
  seq_input{i} = {...
    ['SEQUENCE_ID=' name];
    ['SEQUENCE_TEMPLATE=' seq]; 
    ['SEQUENCE_TARGET=' num2str(targleft(i)) ',' num2str(P.mask_central_region)];
    '=';
  };
end
primer3_input = [primer3_params; cat(1,seq_input{:})];

% switch to temporary directory
tempdir = [P.outstem '.tempdir'];
if ~exist(tempdir,'dir'), mkdir(tempdir); end
p = pwd;
cd(tempdir);

% run primer3
primer3_exe = '/xchip/cga1/lawrence/primer3/primer3-2.2.2-beta/src/primer3_core';
primer3_exe_params = '-strict_tags';
primer3_infile = 'primer3_input.txt';
save_lines(primer3_input,primer3_infile);
primer3_outfile = 'primer3_output.txt';
primer3_errfile = 'primer3_error.txt';
if exist(primer3_outfile,'file'), delete(primer3_outfile); end
if exist(primer3_errfile,'file'), delete(primer3_errfile); end
cmd = [primer3_exe ' ' primer3_exe_params ' -output=' primer3_outfile ' -error=' primer3_errfile ' ' primer3_infile];
fprintf('Running primer3... '); tic
[err res] = system(cmd);
fprintf(' done. '); toc
if err~=0
  disp(load_textfile(primer3_errfile));
  error('primer3 returned nonzero exit code');
end
if ~exist(primer3_outfile,'file'), error('primer3 failed to generate output'); end

% parse output file

p = load_boulder_file(primer3_outfile);

if isfield(p,'PRIMER_ERROR')
  fprintf('PRIMER_ERROR:');
  count(p.PRIMER_ERROR);
  return
end

if isfield(p,'PRIMER_WARNING')
  fprintf('PRIMER_WARNING:');
  count(p.PRIMER_WARNING);
end

flds1 = {'SEQUENCE_ID','SEQUENCE_TEMPLATE','PRIMER_LEFT_0_SEQUENCE','PRIMER_RIGHT_0_SEQUENCE','PRIMER_LEFT_0','PRIMER_RIGHT_0'};
flds2 = {'PRIMER_LEFT_0_TM','PRIMER_RIGHT_0_TM','PRIMER_LEFT_0_GC_PERCENT','PRIMER_RIGHT_0_GC_PERCENT'};
try
  require_fields(p,union(flds1,flds2));
catch me
  fprintf('All designs failed\n');
  return
end

p = make_numeric(p,flds2);
tmp1 = parse(p.PRIMER_LEFT_0,'(\d+),(\d+)',{'leftpos','leftlen'},1:2);
tmp2 = parse(p.PRIMER_RIGHT_0,'(\d+),(\d+)',{'rightpos','rightlen'},1:2);
p = merge_structs({p,tmp1,tmp2});
p.ampstart = p.leftpos;
p.ampend = p.rightpos;
p.amplength = p.ampend - p.ampstart + 1;

% figure(2);draw_primerplot(p);figure(4);hist(p.amplength);

% in-silico PCR

ispcr_input = {};
for i=1:slength(p)
  name = p.SEQUENCE_ID{i};
  f = p.PRIMER_LEFT_0_SEQUENCE{i};
  r = p.PRIMER_RIGHT_0_SEQUENCE{i};
  if isempty(f) || isempty(r), continue; end
  ispcr_input{end+1} = [name ' ' f ' ' r];
end

ispcr_exe = '/broad/tools/Linux/i686/pkgs/isPCR/gfPcr';
ispcr_exe_params = '-maxSize=1000 ispcr 49518 /broad/data/isPCR -out=bed';
ispcr_infile = 'ispcr_input.txt';
save_lines(ispcr_input,ispcr_infile);
ispcr_outfile = 'ispcr_output.txt';
if exist(ispcr_outfile,'file'), delete(ispcr_outfile); end
cmd = [ispcr_exe ' ' ispcr_exe_params ' ' ispcr_infile ' ' ispcr_outfile];

fprintf('Running in-silico PCR... '); tic
[err res] = system(cmd);
fprintf('done. '); toc

if err~=0
  disp(res);
  error('ispcr returned nonzero exit code');
end

if ~exist(ispcr_outfile,'file'), error('isPCR failed to generate output'); end

% parse isPCR output

x = load_lines(ispcr_outfile);
t = parse(x,{'(chr\S*)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)'},{'chr','start','end','id','score','strand'},[2 3 5]);

if slength(t)>0
  fprintf('The following amplicons failed isPCR test\n');
  count(t.id,1);
  keyboard
end


cd(pwd);


% save primers

O = []
O.name = R.seqname;
O.left = p.PRIMER_LEFT_0_SEQUENCE;
O.right = p.PRIMER_RIGHT_0_SEQUENCE;

rfile = [P.outstem '_primers.txt'];
fprintf('Saving %s\n',rfile);
save_struct(O,rfile);
