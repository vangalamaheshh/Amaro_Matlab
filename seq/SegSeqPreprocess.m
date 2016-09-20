function SegSeqPreprocess(sample,P)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 
%   SegSeqPreprocess(sample,P)
%       
%   inputs
%     sample: directory that contains .bam and .bai files
%     old parameter style: SegSeqPreprocess(sample,unpaired_flag)
%
%     P.unpaired_flag
%        if true, tells SegSeqPreprocess.java to accept a BAM with unpaired reads
%         =========================================================================
%           Usage: SegSeqPreprocess <BAMFilename> <BlacklistFilename> <OutFilename> <Chromosome> [mqual] [unpaired]
%                specifying 'mqual' as fifth parameter substitutes mapping qualities for number of mismatches.
%                specifying 'unpaired' as sixth parameter allows to work with unpaired BAMs
%         =========================================================================
%   dependencies
%        trunk/analysis_pipeline/tools/src/java/org/broadinstitute/cga/tools/seq/SegSeqPreprocess.java
%        trunk/analysis_pipeline/genepattern/common_tools/picard/sam-1.03.jar
%   Mike Lawrence 2009, last updated Cheng-Zhong 2011
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('sample','var'), error('<sample> is required'); end

if iscell(sample), error('Multiple samples not yet supported'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'unpaired_flag',false);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %function S = impose_default_value(S,field,value)
% %
% % impose_default_value(S,field,value)
% %
% % S is a struct
% % if field does not exist, it is created with the given default value.
% % if field does exist but is empty, it is given the default value.
% % if field does exist but is nonempty, it is left alone.
% %
% % Mike Lawrence 2008-06-26
% % modified 2009-07-02 to handle required fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


try

java_classpath = [...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar:'...
%   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

short_sample = sample_to_short_sample(sample);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % function short_sample = sample_to_short_sample(sample)
% % short_sample = regexprep(sample,'.*/(.*)/.*','$1');
% % short_sample = regexprep(short_sample,'(.*)-WU','WU-$1');
% % short_sample = regexprep(short_sample,'(.*)-BCM','BCM-$1');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

did_nothing = true;
all_done = false;

while(~all_done)

  all_done = true;

  jobs=[];

  for i=1:2
    if i==1, tn='normal';else tn='tumor'; end
    
    % % Locate paths for inputs and outputs
    
    BAMFile = ['/xchip/tcga_scratch/segseq/' sample '/' tn '.bam'];
    BAIFile = find_bai(BAMFile);

    OutDir = ['/xchip/tcga_scratch/segseq/' sample '/' tn '_SS'];
    
    if ~exist(OutDir,'dir')
      fprintf('Creating directory "%s"\n',OutDir);
      mkdir(OutDir);
    end
    
    dbam = dir(BAMFile);

    if isempty(dbam), error('%s not found',BAMFile); end
    
    % % 
    for c=1:24
      outfile = [OutDir '/chr' num2str(c) '.txt'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout) & dout.bytes>1000
        if dout.datenum >= dbam.datenum, uptodate = true;
        else fprintf('%s need to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'SSP' tn(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' ' ...
           'SegSeqPreprocess ' BAMFile ' ' BAIFile ' ' outfile ' ' num2str(c) ' mqual'];
        if P.unpaired_flag, cmd = [cmd ' unpaired']; end
        cmd = [cmd '"'];
        jobs=[jobs;bsub(cmd,banner)];
      end
    end
    
    outfile = [OutDir '/lanelist.txt'];
    dout = dir(outfile);
    if ~isempty(dout) & dout.bytes>0
      if dout.datenum >= dbam.datenum, uptodate = true;
      else fprintf('%s need to be refreshed:\n',outfile); end
    end
    
    if ~uptodate
      did_nothing = false;
      all_done = false;
      banner = [short_sample 'LANES' tn(1)];
      cmd = ['"java -classpath ' java_classpath ' ' ...
         'MakeLanelist ' BAMFile ' ' outfile '"'];
      jobs=[jobs;bsub(cmd,banner)];
    end
  end % next i   (T/N)

  if ~all_done
    fprintf('Waiting for SegSeqPreprocess to finish\n');
    bwait(jobs);
  end
end % while(~all_done)

if did_nothing
  fprintf('All SegSeq input files already up-to-date.\n');
else
  fprintf('All files now exist.\n');
end

catch me, excuse(me); end






function bainame = find_bai(bamname)


if ~iscell(bamname)
   bainame=test_bam(bamname);
else
  for i=1:length(bamname)
   bainame{i}=test_bam(bamname{i});
  end
end


function bainame = test_bam(bamname)

  input_name=regexprep(bamname,'(.+)\.bam$','$1');
  bamname=[inputname,'.bam'];
  bainame=[inputname,'.bai'];
  if ~exist(bamname), error('bamfile %s does not exist', bamname); end
%  bainame = regexprep(bamname,'\.bam$','\.bai');
%  if strcmp(bainame,bamname) || ~exist(bainame,'file')
%    bainame = [bamname '.bai'];
  if ~exist(bainame,'file'), error('bamfile needs to have an index (bai file) in the same directory with the same name'); end
