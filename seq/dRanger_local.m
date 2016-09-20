function jobs = dRanger_local(sample,P)
% dRanger_local(sample,P)
%

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end
P = impose_default_value(P,'blacklist','none');
P = impose_default_value(P,'window_size',400);
P = impose_default_value(P,'window_step',100);
P = impose_default_value(P,'chr',1:24);

fprintf('dRanger_local\n  sample = %s\n',sample);

java_classpath = [...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

bd = ['/xchip/tcga_scratch/lawrence/' sample];
ss = sample_to_short_sample(sample);

tbam = [bd '/tumor.bam']; demand_file(tbam);
nbam = [bd '/normal.bam']; demand_file(nbam);
tisz = [bd '/tumor.agglom.isz']; demand_file(tisz);
nisz = [bd '/normal.agglom.isz']; demand_file(nisz);
outd = [bd '/dRloc']; if ~exist(outd,'dir'), mkdir(outd); end

jobs = [];
for c=P.chr
  cmd = ['java -Xmx4g -classpath ' java_classpath ' dRangerLocal ' ...
         tbam ' ' nbam ' ' P.blacklist ' ' tisz ' ' nisz ' ' ...
         num2str(P.window_size) ' ' num2str(P.window_step) ' ' ...
         outd '/chr' num2str(c) '.txt ' ...
         num2str(c)];
  banner = [ss 'LOC' num2str(c)];
  jobs = [jobs;bsub(['"' cmd '"'],banner)];
end

