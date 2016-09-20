function X = dRanger_design_primers_for_rearrangements(X,P)
% Mike Lawrence 2010-04-10

if ~exist('P','var'), P = []; end

P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'SNP_mask_BAMs',[]);
P = impose_default_value(P,'maxreads',10000);
P = impose_default_value(P,'Wr',250);                % total window radius
P = impose_default_value(P,'Nr',25);                 % N window radius
P = impose_default_value(P,'mask_central_region',150);

targleft = round((2*P.Wr - P.mask_central_region)/2);
if targleft<1, error('2*P.Wr < P.mask_central_region'); end
mask = [num2str(targleft) ',' num2str(P.mask_central_region)];

if ~isempty(P.SNP_mask_BAMs), require_fields(P.SNP_mask_BAMs,{'individual','bam'}); end

flds1 = {'individual','name'};
flds2 = {'chr1','pos1','str1','chr2','pos2','str2'};
require_fields(X,union(flds1,flds2));
X = make_numeric(X,flds2);
nx = slength(X);
if length(unique(X.name))<nx, error('X.name must be unique'); end

% build validation target sequences

fprintf('Building target sequences... '); tic

if isempty(P.SNP_mask_BAMs), step = 100; else step = 1; end

X.target = cell(nx,1);
for i=1:nx, if ~mod(i,1), fprintf('%d/%d ',i,step); end
  dA = masked_genome_region(X.individual{i},X.chr1(i),X.pos1(i)-P.Wr,X.pos1(i)+P.Wr-1,P.build);
  if X.str1(i)==1, dA = rc(dA); end
  dB = masked_genome_region(X.individual{i},X.chr2(i),X.pos2(i)-P.Wr,X.pos2(i)+P.Wr-1,P.build);
  if X.str2(i)==0, dB = rc(dB); end
  X.target{i} = [dA(1:P.Wr-P.Nr) repmat('N',1,P.Nr*2) dB(P.Wr+P.Nr+1:end)];
end, fprintf('\n');

toc

  function dna = masked_genome_region(individual,chr,st,en,build)
    dna = genome_region(chr,st,en,build);
    bases = 'ACGT';
    if ~isempty(P.SNP_mask_BAMs)
      idx = find(strcmp(individual,P.SNP_mask_BAMs.individual));
      B = cell(length(idx),1);
      if isempty(idx), fprintf('No BAMS match %s\n', individual); return; end
      for bi=1:length(idx)
        [tmp B{bi}] = pull_from_bam(P.SNP_mask_BAMs.bam{idx(bi)},chr,st,en,P);
      end
      B = cat(1,B{~cellfun('isempty',B)});
      if isempty(B), fprintf('BAMS contain no data for this region'); return; end
      range = (st:en)';
      B = B(B(:,2)>25,:);   % keep only high-quality bases
      ntot = as_row(histc(B(:,4),range));
      B = B(B(:,1)>=64,:);  % keep only nonref ACGT
      nmut = as_row(histc(B(:,4),range));
      fracmut = nmut./ntot;
      % number of unique nonref bases at each position
      nnrb = zeros(1,length(range));
      uniquenonrefbase = nan(1,length(range));
      hasdel = false(1,length(range));
      for k=range(1):range(end),j=k-range(1)+1;
        q = unique(B(B(:,4)==k,1));
        nnrb(j) = length(q);
        if length(q)==1, uniquenonrefbase(j) = q; end
        hasdel(j) = ismember(-100,q);
      end
      hznref = (fracmut>=0.9 & ~hasdel & nnrb==1);   % homozygous SNP in tumor+normal
      ismut = find(ntot>=10 & fracmut>=0.15);
      fprintf('\n');
      disp([ismut' ntot(ismut)' nmut(ismut)' nnrb(ismut)' hasdel(ismut)' uniquenonrefbase(ismut)' hznref(ismut)'])
      for k=1:length(ismut),j=ismut(k);
        if hznref(j)
          nonrefbase = bases(uniquenonrefbase(j)-64);
        else
          nonrefbase = 'N';
        end
        dna(j) = nonrefbase;
      end
    end
  end

% PRIMER3 + isPCR

% switch to temporary directory
while(true)
  tempdir = ['/tmp/' num2str(round(mod(now*1e9,1e9)))];
  if ~exist(tempdir,'dir'), mkdir(tempdir); break; end
end
p = pwd;
cd(tempdir);

% design primers

primer3_constant_params = {
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
};

primer3_iteration_params = {
  { 'PRIMER_MIN_GC=40';'PRIMER_MAX_GC=60';'PRIMER_MAX_POLY_X=6' };
  { 'PRIMER_MIN_GC=39';'PRIMER_MAX_GC=60';'PRIMER_MAX_POLY_X=6' };
  { 'PRIMER_MIN_GC=39';'PRIMER_MAX_GC=61';'PRIMER_MAX_POLY_X=6' };
  { 'PRIMER_MIN_GC=38';'PRIMER_MAX_GC=61';'PRIMER_MAX_POLY_X=6' };
  { 'PRIMER_MIN_GC=38';'PRIMER_MAX_GC=62';'PRIMER_MAX_POLY_X=6' };
  { 'PRIMER_MIN_GC=37';'PRIMER_MAX_GC=62';'PRIMER_MAX_POLY_X=7' };
  { 'PRIMER_MIN_GC=37';'PRIMER_MAX_GC=63';'PRIMER_MAX_POLY_X=7' };
  { 'PRIMER_MIN_GC=36';'PRIMER_MAX_GC=63';'PRIMER_MAX_POLY_X=7' };
  { 'PRIMER_MIN_GC=36';'PRIMER_MAX_GC=64';'PRIMER_MAX_POLY_X=7' };
  { 'PRIMER_MIN_GC=35';'PRIMER_MAX_GC=64';'PRIMER_MAX_POLY_X=8' };
  { 'PRIMER_MIN_GC=35';'PRIMER_MAX_GC=65';'PRIMER_MAX_POLY_X=8' };
  { 'PRIMER_MIN_GC=34';'PRIMER_MAX_GC=65';'PRIMER_MAX_POLY_X=8' };
  { 'PRIMER_MIN_GC=34';'PRIMER_MAX_GC=66';'PRIMER_MAX_POLY_X=9' };
  { 'PRIMER_MIN_GC=33';'PRIMER_MAX_GC=66';'PRIMER_MAX_POLY_X=9' };
  { 'PRIMER_MIN_GC=33';'PRIMER_MAX_GC=67';'PRIMER_MAX_POLY_X=9' };
  { 'PRIMER_MIN_GC=32';'PRIMER_MAX_GC=67';'PRIMER_MAX_POLY_X=10' };
  { 'PRIMER_MIN_GC=32';'PRIMER_MAX_GC=68';'PRIMER_MAX_POLY_X=10' };
  { 'PRIMER_MIN_GC=31';'PRIMER_MAX_GC=68';'PRIMER_MAX_POLY_X=10' };
  { 'PRIMER_MIN_GC=31';'PRIMER_MAX_GC=69';'PRIMER_MAX_POLY_X=11' };
  { 'PRIMER_MIN_GC=30';'PRIMER_MAX_GC=69';'PRIMER_MAX_POLY_X=11' };
  { 'PRIMER_MIN_GC=30';'PRIMER_MAX_GC=70';'PRIMER_MAX_POLY_X=11' };
  { 'PRIMER_MIN_GC=30';'PRIMER_MAX_GC=70';'PRIMER_MAX_POLY_X=12' };
  { 'PRIMER_MIN_GC=30';'PRIMER_MAX_GC=70';'PRIMER_MAX_POLY_X=13' };
  { 'PRIMER_MIN_GC=30';'PRIMER_MAX_GC=70';'PRIMER_MAX_POLY_X=14' };
  { 'PRIMER_MIN_GC=30';'PRIMER_MAX_GC=70';'PRIMER_MAX_POLY_X=15' };
};

X.designed = false(nx,1);
X.design_iteration = nan(nx,1);
X.leftprimer = repmat({'none'},nx,1);
X.rightprimer = repmat({'none'},nx,1);
if ~isfield(X,'exclude_regions'), X.exclude_regions = repmat({''},nx,1); end

max_iterations = length(primer3_iteration_params);
for iteration=1:max_iterations
  not_designed = find(~X.designed);

  fprintf('\nIteration %d/%d: designing %d amplicons\n',iteration,max_iterations,length(not_designed));

  seq_input = cell(length(not_designed),1);
  for i=1:length(not_designed),j=not_designed(i);
    seq_input{i} = {
        ['SEQUENCE_ID=' X.name{j}];
        ['SEQUENCE_TEMPLATE=' X.target{j}]; 
        ['SEQUENCE_TARGET=' mask];
        };
    if ~isempty(X.exclude_regions{j})
      seq_input{i} = [seq_input{i};'SEQUENCE_EXCLUDED_REGION=' X.exclude_regions{j}];
    end
    seq_input{i} = [seq_input{i};'='];
  end

  % run primer3

  primer3_params = [primer3_constant_params;primer3_iteration_params{iteration}];
  primer3_input = [primer3_params; cat(1,seq_input{:})];
  primer3_exe = '/xchip/cga1/lawrence/primer3/primer3-2.2.2-beta/src/primer3_core';
  primer3_exe_params = '-strict_tags';
  primer3_infile = 'primer3_input.txt';
  save_lines(primer3_input,primer3_infile);
  primer3_outfile = 'primer3_output.txt';
  primer3_errfile = 'primer3_error.txt';
  if exist(primer3_outfile,'file'), delete(primer3_outfile); end
  if exist(primer3_errfile,'file'), delete(primer3_errfile); end
  cmd = [primer3_exe ' ' primer3_exe_params ' -output=' primer3_outfile ' -error=' primer3_errfile ' ' primer3_infile];
  fprintf('    Running primer3... '); tic
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
    idx = find(~cellfun('isempty',p.PRIMER_ERROR));
    count(p.PRIMER_ERROR(idx));
  end
  
  if isfield(p,'PRIMER_WARNING')
    fprintf('PRIMER_WARNING:');
    idx = find(~cellfun('isempty',p.PRIMER_WARNING));
    count(p.PRIMER_WARNING(idx));
  end
  
  flds1 = {'SEQUENCE_ID','SEQUENCE_TEMPLATE','PRIMER_LEFT_0_SEQUENCE','PRIMER_RIGHT_0_SEQUENCE','PRIMER_LEFT_0','PRIMER_RIGHT_0'};
  flds2 = {'PRIMER_LEFT_0_TM','PRIMER_RIGHT_0_TM','PRIMER_LEFT_0_GC_PERCENT','PRIMER_RIGHT_0_GC_PERCENT'};
  try
    require_fields(p,union(flds1,flds2));
    p = make_numeric(p,flds2);
    tmp1 = parse(p.PRIMER_LEFT_0,'(\d+),(\d+)',{'leftpos','leftlen'},1:2);
    tmp2 = parse(p.PRIMER_RIGHT_0,'(\d+),(\d+)',{'rightpos','rightlen'},1:2);
    p = merge_structs({p,tmp1,tmp2});
    p.ampstart = p.leftpos;
    p.ampend = p.rightpos;
    p.amplength = p.ampend - p.ampstart + 1;
    % figure(2);draw_primerplot(p);figure(4);hist(p.amplength);
    primer3_failed = find(isnan(p.leftpos) | isnan(p.rightpos));
  catch me
    fprintf('    All designs failed primer3.\n');
    continue;   % no need to run isPCR
  end
  
  if length(primer3_failed)>0
    fprintf('    %d/%d amplicons failed in primer3\n',length(primer3_failed),slength(p));
  end

  % in-silico PCR
  
  ispcr_input = {};
  for i=1:slength(p)
    if ismember(i,primer3_failed), continue; end
    name = p.SEQUENCE_ID{i};
    f = p.PRIMER_LEFT_0_SEQUENCE{i};
    r = p.PRIMER_RIGHT_0_SEQUENCE{i};
    ispcr_input{end+1} = [name ' ' f ' ' r];
  end
  
  ispcr_exe = '/broad/tools/Linux/i686/pkgs/isPCR/gfPcr';
  ispcr_exe_params = '-maxSize=1000 ispcr 49518 /broad/data/isPCR -out=bed';
  ispcr_infile = 'ispcr_input.txt';
  save_lines(ispcr_input,ispcr_infile);
  ispcr_outfile = 'ispcr_output.txt';
  if exist(ispcr_outfile,'file'), delete(ispcr_outfile); end
  cmd = [ispcr_exe ' ' ispcr_exe_params ' ' ispcr_infile ' ' ispcr_outfile];
  
  fprintf('    Running in-silico PCR... '); tic
  [err res] = system(cmd);
  fprintf('done. '); toc
  
  if err~=0
    disp(res);
    error('ispcr returned nonzero exit code');
  end
  
  if ~exist(ispcr_outfile,'file'), error('isPCR failed to generate output'); end
  
  % parse isPCR output
  
  x = load_lines(ispcr_outfile);
  I = parse(x,{'(chr\S*)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)'},{'chr','start','end','id','score','strand'},[2 3 5]);
  I.xidx = listmap(I.id,X.name);
  I.pidx = listmap(I.id,X.name(not_designed));
  I.chr = convert_chr(I.chr);
  I.size = I.end-I.start+1;

  % for small indels, some of these isPCR hits might be the actual event
  okidx = find(I.chr==X.chr1(I.xidx) & I.chr==X.chr2(I.xidx) & abs(I.start-X.pos1(I.xidx))<=1000 & abs(I.end-X.pos2(I.xidx))<=1000);
  I = reorder_struct(I,setdiff(1:slength(I),okidx));

  % check which primers failed in isPCR
  bad = unique(I.id);
  fprintf('    %d/%d amplicons failed in isPCR\n',length(bad),slength(p)-length(primer3_failed));

  % check against manually specified forbidden primers
  hit_forbidden_primer = false(slength(p),1);
  if isfield(X,'forbid_leftprimer')
    for i=1:slength(p)
      if strcmpi(p.PRIMER_LEFT_0_SEQUENCE{i},X.forbid_leftprimer(not_designed(i)))
        hit_forbidden_primer(i) = true;
  end,end,end
  if isfield(X,'forbid_rightprimer')
    for i=1:slength(p)
      if strcmpi(p.PRIMER_RIGHT_0_SEQUENCE{i},X.forbid_rightprimer(not_designed(i)))
        hit_forbidden_primer(i) = true;
  end,end,end
  if any(hit_forbidden_primer)
    fprintf('    %d/%d amplicons hit a forbidden primer sequence\n',sum(hit_forbidden_primer),slength(p));
    bad = [bad; p.SEQUENCE_ID(hit_forbidden_primer)];
  end
  
  % for designs that failed isPCR, mark the centers of the primer positions as excluded.
  badno = listmap(bad,p.SEQUENCE_ID);
  for j=1:length(badno)
    pidx = badno(j); xidx = not_designed(pidx);
    lpos = round(p.leftpos(pidx)+p.leftlen(pidx)/2); rpos = round(p.rightpos(pidx)-p.rightlen(pidx)/2);
    excl = [num2str(lpos) ',1 ' num2str(rpos) ',1 '];
    if isempty(X.exclude_regions{xidx}) X.exclude_regions{xidx} = excl;
    else X.exclude_regions{xidx} = [X.exclude_regions{xidx} excl]; end
  end

  % designs that passed both primer3 and isPCR
  idx = setdiff(1:slength(p),[primer3_failed;I.pidx;find(hit_forbidden_primer)]);
  if ~isempty(idx)
    didx = not_designed(idx);
    X.designed(didx) = true;
    X.design_iteration(didx) = iteration;
    X.leftprimer(didx) = p.PRIMER_LEFT_0_SEQUENCE(idx);
    X.rightprimer(didx) = p.PRIMER_RIGHT_0_SEQUENCE(idx);
    X.amplicon_start(didx) = p.ampstart(idx);
    X.amplicon_end(didx) = p.ampend(idx);
    X.amplicon_length(didx) = p.amplen(idx);
    X.amplicon_sequence(didx) = X.target{didx}(p.ampstart(idx):p.ampend(idx)); 
  end

  if all(X.designed), break; end

end % next iteration

if all(X.designed)
  fprintf('\nAll designed succeeded\n');
else
  fprintf('\n%d/%d designs failed\n',sum(~X.designed),nx);
end

cd(pwd);


end % main function

