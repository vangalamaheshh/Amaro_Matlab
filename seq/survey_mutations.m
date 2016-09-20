function survey_mutations(mutlist,bamfile,blacklist,refdir,outname,P);
%  Mike Lawrence 2009-11-16

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Last parameter should be a P struct'); end

P = impose_default_value(P,'use_LSF',true);
P = impose_default_value(P,'banner',[]);

java_classpath = [...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
%   '/seq/software/picard/current/bin/sam-1.05.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

if ~isempty(P.banner), banner = P.banner; else banner = ''; end
cmd = ['java -classpath ' java_classpath ' ' ...
       'SurveyMutations ' mutlist ' ' bamfile ' ' blacklist ' ' refdir ' ' outname];
if P.use_LSF
  bsub(['"' cmd '"'],banner);
else
  system(cmd);
end
