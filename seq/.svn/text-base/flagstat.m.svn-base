function flagstat(samples,P)
% Mike Lawrence 2009

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);


try

if P.cancer_sample, tn = {'normal','tumor'}; else tn = {'sample'}; end

did_nothing = true;
all_done = false;
while(~all_done)
  jobs=[]; all_done = true;
  for i=1:length(samples)
    sample = samples{i};
    direc = ['/xchip/tcga_scratch/lawrence/' sample];
    for t=1:length(tn)
      bam = [direc '/' tn{i} '.bam'];
      dbam = dir(bam);
      if isempty(dbam)
        fprintf('ERROR: %s not found!\n',bam);
        continue
      end
      out = [direc '/' tn{i} '.stats'];
      dout = dir(out);
      if ~isempty(dout)
        %fprintf('out=%s outdate=%s bamdate=%s\n',out,dout.date,dbam.date);
        if dout.datenum >= dbam.datenum
          continue;   % stats already up-to-date
      end,end
      outtmp = [out '_tmp'];
      all_done = false;
      did_nothing = false;
      cmd = ['"samtools flagstat ' bam ' > ' outtmp ';'...
               'mv ' outtmp ' ' out '"'];
      banner = [sample_to_short_sample(sample) 'FLAG' tn{i}(1)];
      jobs = [jobs; bsub(cmd,banner)];
  end,end
  if (~all_done), bwait(jobs); end
end

if did_nothing, fprintf('Stats files already up-to-date for all those samples.\n'); end

catch me, excuse(me); end
