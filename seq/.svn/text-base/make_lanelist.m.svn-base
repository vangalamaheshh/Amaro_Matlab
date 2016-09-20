function make_lanelist(samples,P)
% (1) deletes any lanelists that are older than the corresponding BAM file
% (2) runs MakeLanelist.java for any BAMS without a lanelist

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'use_LSF',true);

java_classpath = [...
   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

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
      stem = [direc '/' tn{t}];
      bam = [stem '.bam'];
      bai = find_bai(bam);
      lanelist = [stem '.bam.lanelist'];
      dbam = dir(bam);
      if isempty(dbam)
        fprintf('ERROR: %s not found!\n',bam);
        continue
      end
      dlanelist = dir(lanelist);
      if ~isempty(dlanelist)
        if datenum(dlanelist.date)>=datenum(dbam.date)
          continue;          % lanelist is already up-to-date
      end,end
      all_done = false;
      did_nothing = false;
      cmd = ['java -classpath ' java_classpath ' MakeLanelist ' bam ' ' lanelist];
      banner = [sample_to_short_sample(sample) 'LANES' tn{t}(1)];
      if P.use_LSF
        jobs = [jobs; bsub(['"' cmd '"'],banner)];
      else
        system(cmd);
      end
  end,end
  if (~all_done), bwait(jobs); else break; end
end

if did_nothing, fprintf('Lanelists are already up-to-date for all those samples.\n'); end

catch me, excuse(me); end

