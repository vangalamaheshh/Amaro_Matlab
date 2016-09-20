function segseq_run(pat,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'use_LSF',false);
P = impose_default_value(P,'bsub_queue','cga');
P = impose_default_value(P,'queue',P.bsub_queue);
P = impose_default_value(P,'bsub_mem',16);
%P = impose_default_value(P,'csh_location','/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/SegSeqRun.csh');
P = impose_default_value(P,'exe_dir','/xchip/cga2/lawrence/cga/trunk/segseq/preprocess3');
P = impose_default_value(P,'offset_mismatch',0);
P = impose_default_value(P,'qual_cutoff',20);
P = impose_default_value(P,'max_chr',23);
P = impose_default_value(P,'output_all_flowcells',1);

require_fields(pat,{'name','tumor_ss','normal_ss','segseq_basedir'});

cmds = {};
banners = {};
jobs = [];
for i=1:slength(pat)
  if P.use_LSF
    cmd = ['matlab -r "segseq_run_sample(''' pat.normal_ss{i} ''',''' pat.tumor_ss{i} ''',''' pat.name{i} ''',''' ...
         num2str(P.offset_mismatch) ''',''' num2str(P.qual_cutoff) ''',''' P.exe_dir ''',''' ...
         pat.segseq_basedir{i} ''',''' num2str(P.max_chr) ''',''' num2str(P.output_all_flowcells) ''')"'];
    cmds{end+1} = ['-R "rusage[matlab=' num2str(P.bsub_mem) ']" ' cmd];
    banners{end+1} = [pat.name{i} '_SS'];
    % submit jobs as they accumulate
    if length(cmds)>=60
      jobs = [jobs; bsub(cmds,banners,P)];
      cmds = {}; banners = {};
    end
  else
    segseq_run_sample(pat.normal_ss{i},pat.tumor_ss{i},pat.name{i},P.offset_mismatch,P.qual_cutoff,P.exe_dir,...
                      pat.segseq_basedir{i},num2str(P.max_chr),num2str(P.output_all_flowcells));
  end
end

if ~isempty(cmds)
  jobs = [jobs; bsub(cmds,banners,P)];
  cmds = {}; banners = {};
end
    
if P.use_LSF
  fprintf('Waiting for jobs to finish\n');
  bwait(jobs);
end
