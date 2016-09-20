function [C J S] = get_calls_at_positions(bamfiles,chr,pos,outdir,blacklist,qualcutoff,P)
% get_calls_at_positions(bamfiles,chr,pos,outdir,blacklist,qualcutoff)
%
% C = calls (rows=positions, columns=[A C G T DEL INS], pages=bamfiles)
% J = judgement (rows=positions, columns=bamfiles)
%      0 = reference
%     +1 = has non-reference base
%     +2 = has deletion
%     +4 = has insertion
% S = summary (rows=positions, columns=#bams showing evidence of [SNP DEL INS]

if nargin==4 && isstruct(outdir)
  P=outdir;
  outdir=[];
  blacklist = 'none';
  qualcutoff = 0;
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'banner','GCP');
P = impose_default_value(P,'force_report',false);
P = impose_default_value(P,'load_all_results_into_memory',true);
P = impose_default_value(P,'write_binary_file',false);
if isfield(P,'suppress_bamidx') && ~isfield(P,'suppress_bamnames'), error(' please supply P.suppress_bamnames'); end
P = impose_default_value(P,'suppress_bamnames',{});
if isempty(P.suppress_bamnames), P.suppress_bamnames = {}; end % because ismember will choke on []
P = impose_default_value(P,'include_duplicate_reads',false);
if ischar(P.include_duplicate_reads)
  if strcmpi(P.include_duplicate_reads,'true')
    include_dups = '1';
  elseif  strcmpi(P.include_duplicate_reads,'false')
    include_dups = '0';
  else
    error('P.include_duplicate_reads should be true or false');
  end
elseif isnumeric(P.include_duplicate_reads)
  if P.include_duplicate_reads==1 
    include_dups = '1';
  elseif P.include_duplicate_reads==0
    include_dups = '0';
  else
    error('P.include_duplicate_reads should be 1 or 0');
  end
elseif islogical(P.include_duplicate_reads)
  include_dups = num2str(P.include_duplicate_reads);
else
  error('invalid P.include_duplicate_reads');
end

if nargin==5 && ischar(blacklist) && strncmp(blacklist,'hg',2)
  % get_calls_at_positions(bamfiles,chr,pos,outdir,build)
  P.build = blacklist;
  clear blacklist;
end

if ~P.force_report
  if isfield(P,'refdir')
    refdir = P.refdir;
  else
    if ~isfield(P,'build')
      fprintf('Assuming build hg19\n');
%      fprintf('Is this OK? (dbcont/dbquit)\n'); keyboard
      buildno = 19;
    else
      buildno = interpret_build(P.build);
    end
    if buildno==18 || buildno==19
      fprintf(['using default hg' num2str(buildno) ' refdir\n']);
      refdir = ['/xchip/cga1/annotation/db/ucsc/hg' num2str(buildno)];
    else
      error('unknown build %s',P.build);
    end
  end
end

if nargin<3, error('needs at least 3 input args'); end

if nargin==3 || isempty(outdir)   % no outdir assigned
  outdir = ['/xchip/cga1/temp/getcalls' num2str(round(1e12*rand))];
  fprintf('Outputting to temp directory:  %s\n',outdir);
end
ede(outdir);

if ischar(bamfiles), bamfiles = {bamfiles}; end
if ~iscell(bamfiles), error('bamfiles should be a list of BAM files'); end

chr = as_column(chr);
pos = as_column(pos);
if length(chr)==1 && length(pos)>1, chr = chr*ones(size(pos)); end
if length(chr)~=length(pos), error('length(chr)~=length(pos)'); end
nmut = length(chr);

if iscell(chr), chr = convert_chr(chr); end
if ~isnumeric(pos), pos = str2double(pos); end

bad_idx = find(isnan(chr) | isnan(pos));
if ~isempty(bad_idx)
  fprintf('%d/%d positions had NaN for chr and/or pos:  ignoring these positions\n',length(bad_idx),nmut);
end
good_idx = setdiff(1:nmut,bad_idx);

if ~exist('blacklist','var'), blacklist = 'none'; end
if ~exist('qualcutoff','var'), qualcutoff = 20; end

loaded_saved_results = false;
finished_file = [outdir '/C.mat'];
if exist(finished_file,'file')
  fprintf('Loading finished results from %s\n',finished_file);
  load(finished_file,'C','ref_base');
  loaded_saved_results = (exist('C','var') & exist('ref_base','var'));
end

if ~loaded_saved_results
  if ~P.force_report
    % for output files that already exist, don't worry about missing BAM/BAI's
    for f=1:length(bamfiles)
      outfile_already_exists(f,1) = ~~exist([outdir '/calls.' num2str(f) '.txt'],'file');
    end

    if any(outfile_already_exists)
      fprintf('%d/%d output files already exist\n',sum(outfile_already_exists),length(bamfiles));
    end

    if exist([outdir '/bamfile_list.txt'],'file')
      % bamfile list is already finalized: do not alter
      bamfiles = load_lines([outdir '/bamfile_list.txt']);
    else
      % remove non-existent files
      btmp = bamfiles(~outfile_already_exists);
      ok = demand_files(btmp);
      if any(~ok)
        fprintf('WARNING: %d/%d files do not exist, will be skipped.\n',sum(~ok),length(bamfiles));
        bamfiles(ismember(bamfiles,btmp(~ok))) = [];
      end

      % make sure BAM indices are also available
      baifiles = regexprep(bamfiles,'\.bam$','.bai');
      itmp = baifiles(~outfile_already_exists);
      btmp = bamfiles(~outfile_already_exists);
      ok = demand_file(itmp);
      if any(~ok)
        itmp2 = itmp(~ok);
        btmp2 = btmp(~ok);
        itmp2 = regexprep(itmp2,'\.bai$','.bam.bai');
        ok = demand_files(itmp2);
        if any(~ok)
          fprintf('WARNING: %d/%d index files do not exist; those BAMs will be skipped.\n',sum(~ok),length(bamfiles));
          bamfiles(ismember(bamfiles,btmp2(~ok))) = [];
        else
          fprintf('But *.bam.bai files ARE available in all of the above cases.\n');
        end
      end
    end

    save_lines(bamfiles,[outdir '/bamfile_list.txt']);
  end

  for f=1:length(bamfiles)
    outfiles{f,1} = [outdir '/calls.' num2str(f) '.txt'];
    tmpfiles{f,1} = [outdir '/calls.' num2str(f) '.txt_partial'];
  end
  
  X = [];
  X.chr = chr(good_idx);
  X.pos = pos(good_idx);
  poslist = [outdir '/position_list.txt'];
  if ~exist('poslist','file'), save_struct_noheader(X,poslist); end

  if ~P.force_report
    % SUBMIT JOBS
    java_classpath = get_jcp();
    java_classname = 'org/broadinstitute/cga/tools/seq/GetCallsAtPositions';
    
    did_nothing = true;
    all_done = false;
    while(~all_done)
      all_done = true;
      jobs=[];
      cmds = {}; banners = {};
      for f=1:length(bamfiles)
        if ~exist(outfiles{f},'file') && ~ismember(bamfiles{f},P.suppress_bamnames) %&& f<2900
          did_nothing = false;
          all_done = false;
          banners{end+1,1} = [P.banner num2str(f)];
          cmds{end+1} = ['"java -classpath ' java_classpath ' ' java_classname ' ' ...
                         bamfiles{f} ' ' poslist ' ' outfiles{f} ' 3 ' num2str(qualcutoff) ...
                         ' ' blacklist ' ' refdir ' ' include_dups '"'];
          % submit jobs as they accumulate
          if length(cmds)>=60
            jobs = [jobs;bsub(cmds,banners,P)];
            cmds = {}; banners = {};
          end
        end
      end
      if ~isempty(cmds)
        jobs = [jobs;bsub(cmds,banners,P)];
        cmds = {}; banners = {};
      end
    
      if ~all_done
        fprintf('Waiting for jobs to finish\n');
        bwait(jobs);
      end
    end % while(~all_done)
    
    ok = demand_files(outfiles);
    
    if did_nothing
      fprintf('All files already up-to-date.\n');
    else
      fprintf('All files now exist.\n');
    end
    
  end  % if ~P.force_report

  % load results
  fprintf('Loading results: ');

  for f=1:length(outfiles)
    try
      ncols = get_colcount(outfiles{f},0);
      break
    catch me
    end
  end
  if ~exist('ncols','var'), error('No output files exist'); end
  if ncols~=7 && ncols~=9, error('wrong column count in first result file'); end
  fmt = ['%f%f%s' repmat('%f',1,ncols-3)];
  
  if P.load_all_results_into_memory
    try 
      D = zeros(slength(X),ncols-3,length(outfiles));
    catch me
      fprintf('Not enough memory to use "double"... trying "int16"\n');
      D = zeros(slength(X),ncols-3,length(outfiles),'int16');
    end 
  end
  
  if P.write_binary_file
    binfilename = [outdir '/results.bin'];
    fprintf('Writing binary file %s\n',binfilename);
    BIN = fopen(binfilename,'w');
    block = zeros(nmut*6,1,'uint8');
    fprintf('Blocksize = %d\n', nmut*6);
  end

  clear ref_base
  for f=1:length(outfiles), if ~mod(f,10), fprintf('%d/%d ',f,length(outfiles)); end
    try
      if exist(outfiles{f},'file') && ~ismember(bamfiles{f},P.suppress_bamnames)
        tmp = read_table(outfiles{f},fmt,char(9),0);
        if ~exist('ref_base','var'), ref_base = tmp.dat{3}; end
        rawblock = cat(2,tmp.dat{4:end});
        if size(rawblock,1)~=slength(X), error('Wrong number of mutations'); end
        if size(rawblock,2)~=ncols-3, error('Wrong number of columns'); end
        if P.load_all_results_into_memory
          D(:,:,f) = rawblock;
        end
        if P.write_binary_file
          block(:) = rawblock(:);
          block(rawblock>255) = 255;
        end
      else
        fprintf('Skipping nonexistent output file %d\n',f);
      end
    catch me
      if P.write_binary_file
        block(:)=0;
      end
      fprintf('Trouble loading file %d: ',f);
      disp(me.message);
%      keyboard
    end
    if P.write_binary_file
      fwrite(BIN,block,'uint8');
    end
  end, fprintf('\n');
  
  if P.write_binary_file
    fclose(BIN);
  end

  if ~P.load_all_results_into_memory
    fprintf('Wrote binary file, did not load all results into memory.\n');
    C=[];
    J=[];
    S=[];
    return;
  end
  
  % restore "bad" positions
  try
    C = zeros(nmut,ncols-3,length(outfiles));
  catch me
    fprintf('Not enough memory to use "double"... trying "int16"\n');
    C = zeros(nmut,ncols-3,length(outfiles),'int16');
  end
  C(good_idx,:,:) = D;

  try
    save(finished_file,'C','ref_base');
  catch me
    disp(me)
    disp(me.message)
    fprintf('Failed to save %s, probably because C too large\n', finished_file);
  end

end

% summarize results for convenience

C = C(good_idx, :, :);

if nargout>1
  nr = sum(C,2);
  C2 = bsxfun(@rdivide,C,nr);

  refi = listmap(ref_base,{'A','C','G','T'});
  if sum(isnan(refi))>=0.5 && nargout>1
    fprintf('WARNING:  Reference data needed to determine SNP status in returned matrix J\n');
  end

  C3 = C2(:,1:4,:); for i=1:size(C3,1), if ~isnan(refi(i)), C3(i,refi(i),:)=0; end; end
  nonref = sum(C3,2);
  
  has_snp = squeeze((nr>=20 & nonref>=0.2));
  has_del = squeeze((nr>=10 & C2(:,5,:)>=0.1));
  has_ins = squeeze((nr>=10 & C2(:,6,:)>=0.1));
  
  J = has_snp + 2*has_del + 4*has_ins;

  % display as heatmap
  [tmp row_ord] = sort(sum(J,2));
  [tmp col_ord] = sort(sum(J,1));
  imagesc(J(row_ord,col_ord));
  
  nsnp = sum(has_snp,2);
  ndel = sum(has_del,2);
  nins = sum(has_ins,2);
  S = [nsnp ndel nins];

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% further stuff you can do if you know the type of mutation you're looking for
% --a lot of this is now done in survey_panel_of_normals_for_mutations.m

% overall judgement per mutation
M.survey = S;
M.is_ins = grepmi('ins',M.type);
M.is_del = grepmi('del',M.type);
M.is_bad_ins = (M.is_ins & S(:,3)>=2);     % 2 samples enough to throw out an insertion
M.is_bad_del = (M.is_del & S(:,2)>=2);     % 2 samples enough to throw out a deletion
M.is_bad_snp = (M.is_point & S(:,1)>=5);   % 5 samples needed to throw out a point mutation
M.is_bad = (M.is_bad_ins | M.is_bad_del | M.is_bad_snp);

