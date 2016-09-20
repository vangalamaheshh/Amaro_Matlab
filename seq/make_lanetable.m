function make_lanetable(samples,P)
% (1) deletes any lanelists that are older than the corresponding BAM file
% (2) runs MakeLanetable.java for any BAMS without a lanetable
% (3) then for each file, adds a "baitset" column


if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'picard_dir','/seq/picard/by_project/C239');
P = impose_default_value(P,'use_LSF',true);

java_classpath = [...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
%   '/seq/software/picard/current/bin/sam-1.05.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

try

if P.cancer_sample, tn = {'normal','tumor'}; else tn = {'sample'}; end

did_nothing = true;
all_done = false;
files_updated = {};
while(~all_done)
  jobs=[]; all_done = true;
  for i=1:length(samples)
    sample = samples{i};
    direc = ['/xchip/tcga_scratch/lawrence/' sample];
    for t=1:length(tn)
      stem = [direc '/' tn{t}];
      bam = [stem '.bam'];
      lanetable = [stem '.bam.lanetable'];
      dbam = dir(bam);
      if isempty(dbam)
        fprintf('ERROR: %s not found!\n',bam);
        continue
      end
      dlanetable = dir(lanetable);
      if ~isempty(dlanetable)
        if datenum(dlanetable.date)>=datenum(dbam.date)
          continue;          % lanetable is already up-to-date
      end,end
      all_done = false;
      did_nothing = false;
      files_updated = [files_updated; lanetable];
      cmd = ['java -classpath ' java_classpath ' MakeLanetable ' bam ' ' lanetable];
      banner = [sample_to_short_sample(sample) 'LNTAB' tn{t}(1)];
      if P.use_LSF
        jobs = [jobs; bsub(['"' cmd '"'],banner)];
      else
        system(cmd);
      end
  end,end
  if (~all_done), bwait(jobs); else break; end
end

if did_nothing
  fprintf('Lanetables are already up-to-date for all those samples.\n');
  return
end

% add "baitset" column

fprintf('Loading library <-> baitset correspondences\n');
X = get_library_list(P.picard_dir);
fprintf('Annotating baitsets in file: ');
for i=1:length(files_updated),fprintf('%d/%d ',i,length(files_updated));
  fname = files_updated{i};
  T = load_struct(fname);
  T.baitset = mapacross(T.LB,X.library,X.baitset,'unknown');
  save_struct(T,fname);
end, fprintf('\n');

catch me, excuse(me); end
