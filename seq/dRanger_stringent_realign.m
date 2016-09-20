function dRanger_stringent_realign(sample,P)
% dRanger_stringent_realign(sample,P)
%
% Mike Lawrence 2009-07-23

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Multiple samples not supported'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'tumor_only',false);
P = impose_default_value(P,'normal_only',false);
P = impose_default_value(P,'vicinity_forward',10000);
P = impose_default_value(P,'vicinity_backward',100);
P = impose_default_value(P,'self_align_margin',20);
P = impose_default_value(P,'max_tasks_per_batch',50000);
P = impose_default_value(P,'pause_before_submitting_batches',false);
P = impose_default_value(P,'pause_after_submitting_batches',false);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
P = impose_default_value(P,'dRanger_input_matfile','dRanger_input.mat');
P = impose_default_value(P,'dRanger_seqs_matfile','dRanger_seqs.mat');
P = impose_default_value(P,'working_subdir','sa_4');
P = impose_default_value(P,'java_classname','BatchAlign2');
P = impose_default_value(P,'scoring_matrix','/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/NT3.txt');
P = impose_default_value(P,'gap_open_penalty',5);
P = impose_default_value(P,'gap_extend_penalty',1);


%%
%%
fprintf('2009-10-15\n');
fprintf('dRanger_stringent_realign needs to be modified before use:\n');
fprintf('Input and output files need to be changed\n');
keyboard;
%%
%%


java_classpath = [...
   '/xchip/tcga_scratch/lawrence/jaligner/jaligner-1.0/jaligner.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

short_sample = sample_to_short_sample(sample);

len = load_chrlen('chrfiles');

if P.cancer_sample, tn = {'normal','tumor'}; else tn = {'sample'}; end

T = cell(length(tn),1);
B = cell(length(tn),1);

try

did_nothing = true;
all_done = false;
while(~all_done)
 all_done = true;
 jobs = [];

 for i=1:length(tn)
  if strcmp(tn{i},'normal') & P.tumor_only, continue; end
  if strcmp(tn{i},'tumor') & P.normal_only, continue; end

  fprintf('PROCESSING %s\n',upper(tn{i}));
  direc = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '_' P.dRangerPreprocess_output_dir_suffix];

  % see if output file exists and doesn't need refreshing
  infile = [direc '/' P.dRangerPreprocess_pairs_file];
  din = dir(infile);
  if isempty(din), error('Not found: %s',infile); end
  outfile = [direc '/' P.dRanger_stringent_realign_output_file];
  dout = dir(outfile);
  if ~isempty(dout)
    fprintf('stringent_realign output file already exists ');
    if dout.datenum>=din.datenum
      fprintf('and does not need to be refreshed\n');
      continue;
    else
      fprintf('but needs to be refreshed\n');
    end
  end

  sadirec = [direc '/' P.working_subdir];
  if ~exist(sadirec,'dir'), mkdir(sadirec); end

  if isempty(T{i})
    % load input file
    fprintf('Loading %s\n',infile);
    tic, X = load_struct(infile,[repmat('%f',1,13) '%s%s'],0); toc
    X = rename_fields(X,colx(1:15),{'rgrp','namenumber','chr1','start1','end1','strand1','qual1',...
                      'chr2','start2','end2','strand2','qual2','flip','seq1','seq2'});
    nx = slength(X);

    % reformat into tasks:  end1->self end2->self end1->other end2->other
    fprintf('Reformatting into tasks\n');
    T{i}.record = repmat((1:nx)',4,1);
    T{i}.whichend = repmat([ones(nx,1);2*ones(nx,1)],2,1);
    T{i}.self = [ones(2*nx,1); zeros(2*nx,1)];
    T{i}.seq = repmat([X.seq1; X.seq2],2,1);
    T{i}.rcflag = repmat([X.strand1; X.strand2],2,1);
    T{i}.vchr = [X.chr1; X.chr2; X.chr2; X.chr1];
    ostrand = 1 - [X.strand2; X.strand1];
    ostart = [X.start2; X.start1] - P.vicinity_backward;
    tmp = [X.end2; X.end1] - P.vicinity_forward;
    ostart(ostrand==0) = tmp(ostrand==0);
    oend = ostart + P.vicinity_forward + P.vicinity_backward;
    T{i}.vstrand = [[X.strand1; X.strand2]; ostrand];
    T{i}.vstart = [[X.start1; X.start2] - P.self_align_margin; ostart];
    T{i}.vend = [[X.end1; X.end2] + P.self_align_margin; oend];

    % clear X nx;

    % remove tasks that we don't have sequence info for (i.e. at ends of chromosomes)
    T{i}.vstart = max(1,T{i}.vstart);
    T{i}.vstart = min(len(T{i}.vchr),T{i}.vstart);
    T{i}.vend = min(len(T{i}.vchr),T{i}.vend);
    T{i} = reorder_struct(T{i},T{i}.vend>T{i}.vstart);
    nt = slength(T{i});

    % sort by vicinity: chromosome, strand, start, end
    fprintf('Sorting tasks\n');
    T{i} = sort_struct(T{i}, {'vchr','vstrand','vstart','vend'});
  end

  if isempty(B{i})
    % divide tasks into batches
    fprintf('Dividing tasks into batches\n');
    B{i}.first_task = []; B{i}.last_task = [];
    flip_Tvchr = flipud(T{i}.vchr);
    flip_Tvstrand = flipud(T{i}.vstrand);
    step = P.max_tasks_per_batch;
    for vchr=1:24
      for vstrand=0:1
        f = find(T{i}.vchr==vchr & T{i}.vstrand==vstrand, 1);
        if ~isempty(f)
          l = nt+1-find(flip_Tvchr==vchr & flip_Tvstrand==vstrand, 1);
          ff = f:step:l-1;
          ll = ff+(step-1); ll(end) = l;
          B{i}.first_task = [B{i}.first_task; ff'];
          B{i}.last_task = [B{i}.last_task; ll'];
        end
      end
    end
    B{i}.vchr = T{i}.vchr(B{i}.first_task);
    B{i}.vstrand = T{i}.vstrand(B{i}.first_task);
    B{i}.vmin = T{i}.vstart(B{i}.first_task);
    B{i}.vmax = T{i}.vend(B{i}.last_task);
    nb(i) = slength(B{i});
    bfile = [sadirec '/batchlist.txt'];
    tmp = []; tmp.num = (1:nb(i))';
    tmp = merge_structs({tmp,B{i}});
    save_struct(tmp,bfile);
  end

  if P.pause_before_submitting_batches
    fprintf('Type "return" to begin submitting batches.\n');
    keyboard
  end

  % ALIGN: write inputfiles for Java aligner and dispatch batches

  for b=1:nb(i)
    if B{i}.vstrand(b)==0, sc = '+'; else sc = '-'; end
    fprintf('%s batch %d/%d chr%d:%d-%d(%s)\n',upper(tn{i}),b,nb(i),...
       B{i}.vchr(b),B{i}.vmin(b),B{i}.vmax(b),sc);

    % vicinity file
    vfile = [sadirec '/vicinity_batch' num2str(b) '.txt'];    dirv = dir(vfile);
    uptodate = false;
    if ~isempty(dirv)
      if dirv.datenum>=din.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',vfile); end
    end
    if ~uptodate
      fprintf('Writing vicinity DNA to %s\n',vfile);
      vdna = genome_region(B{i}.vchr(b),B{i}.vmin(b),B{i}.vmax(b));
      if B{i}.vstrand(b)==1, vdna = rc(vdna); end
      save_textfile(vdna,vfile);
    end

    % sequences file
    sfile = [sadirec '/seqs_batch' num2str(b) '.txt'];
    dirs = dir(sfile);
    uptodate = false;
    if ~isempty(dirs)
      if dirs.datenum>=din.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',sfile); end
    end
    if ~uptodate
      fprintf('Writing sequences DNA to %s\n',sfile);
      S = [];
      tasks = B{i}.first_task(b):B{i}.last_task(b);
      S.record = T{i}.record(tasks);
      S.whichend = T{i}.whichend(tasks);
      S.self = T{i}.self(tasks);
      S.seq = T{i}.seq(tasks);
      S.rcflag = T{i}.rcflag(tasks);
      S.startrel = T{i}.vstart(tasks) - B{i}.vmin(b);
      S.endrel = T{i}.vend(tasks) - B{i}.vmin(b);
      if B{i}.vstrand(b)==1
        sz = B{i}.vmax(b) - B{i}.vmin(b) + 1;
        tmp = S.startrel;
        S.startrel = sz - S.endrel;
        S.endrel = sz - tmp;
      end
      save_struct(S,sfile,'no_headers');
    end
   
    % call BatchAlign
    rfile = [sadirec '/results_batch' num2str(b) '.txt'];
    dirr = dir(rfile);
    uptodate = false;
    if ~isempty(dirr)
      if dirr.datenum>=din.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',rfile); end
    end
    if ~uptodate
      all_done = false;
      did_nothing = false;
      banner = [short_sample 'BA' tn{i}(1) num2str(b)];
      cmd = ['"java -Xmx2g -classpath ' java_classpath ' ' P.java_classname ' ' ...
             P.scoring_matrix ' '  vfile ' ' sfile ' ' rfile ' ' ...
             num2str(P.gap_open_penalty) ' ' num2str(P.gap_extend_penalty) ...
             ' >/dev/null 2>&1"'];
      jobs = [jobs;bsub(cmd,banner)];
    end
  end    % next b  (batch)

  if P.pause_after_submitting_batches
    fprintf('Finished with %s.  Type "return" to continue.\n',upper(tn{i}));
    keyboard
  end

 end % next i  (T/N)

 if ~all_done
   fprintf('Waiting for BatchAlign jobs to finish\n');
   bwait(jobs);
 else
   break
 end

end % while(~all_done)

if did_nothing
  fprintf('All dRanger_stringent_realign results files are already up-to-date.\n');
else
  fprintf('All dRanger_stringent_realign results files are now up-to-date.\n');
end

dRanger_stringent_realign_process_results(sample,P);

catch me; excuse(me); end
