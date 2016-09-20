function dRanger_bridgepairs(sample,window_size)
if ~exist('window_size','var'),window_size=3000;end

try

fprintf('dRanger_bridgepairs\n\tsample = %s\n\n',sample);

fprintf('Loading dRanger data\n');
direc = ['/xchip/tcga_scratch/lawrence/' sample];
name2 = upper(regexprep(sample,'/','-'));
fname = [direc '/' name2 '_dRanger_results_bfiltered.txt'];
if ~exist(fname,'file'), error('Can''t find file %s',fname);end
X = load_struct(fname);
X = reorder_struct(X,strcmp(X.filterB,'0'));
X = reorder_struct(X,strcmp(X.filterHCL,'0'));
X = make_numeric(X,{'chr1','chr2','pos1','pos2'});

% working directory

workdirT = [direc '/tumor.bridge'];
workdirN = [direc '/normal.bridge'];
if ~exist(workdirT), mkdir(workdirT); end
if ~exist(workdirN), mkdir(workdirN); end

% make target list

T = [];
T.num = (1:slength(X)*2)';
T.chr = [X.chr1;X.chr2];
T.st = [X.pos1;X.pos2]-round(window_size/2);
T.en = [X.pos1;X.pos2]+round(window_size/2);
nt = slength(T);

targfileT = [workdirT '/targs.txt'];
save_struct(T,targfileT,'no_headers');
targfileN = [workdirN '/targs.txt'];
save_struct(T,targfileN,'no_headers');

% call ExtractRegions

cmdT = ['java -classpath /xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
    '/xchip/tcga/gbm/analysis/lawrence/sam ExtractRegions '...
    direc '/tumor.bam ' direc '/tumor.bam.bai ' targfileT ' ' workdirT '/reg'];
cmdN = ['java -classpath /xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
    '/xchip/tcga/gbm/analysis/lawrence/sam ExtractRegions '...
    direc '/normal.bam ' direc '/normal.bam.bai ' targfileN ' ' workdirN '/reg'];
jobs = [bsub(['"' cmdT '"']);bsub(['"' cmdN '"'])];
bwait(jobs);

% sort and join

jobs=[];
for i=1:nt
  cmdT = ['/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sort_and_join '...
    workdirT '/reg' num2str(i) '.txt ' workdirT '/reg' num2str(i) '.sj.txt'];
  cmdN = ['/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sort_and_join '...
    workdirN '/reg' num2str(i) '.txt ' workdirN '/reg' num2str(i) '.sj.txt'];
  jobs = [jobs;bsub(['"' cmdT '"']);bsub(['"' cmdN '"'])];
end
bwait(jobs);



catch me, excuse(me); end
