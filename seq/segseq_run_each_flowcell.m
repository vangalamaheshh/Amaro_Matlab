function segseq_run_each_flowcell(base_dir,fcidx)

dr = [base_dir '/matfiles'];
cd(dr);

d = dir('*fc_tumor_normal_aligned_paired_reads*mat');
flowcellfile = list2cell(d.name);
flowcell = regexprep(flowcellfile,'(.*)_(median|nonmed)fc_tumor_normal_aligned_paired_reads.*','$1');

if ~exist('fcidx','var'), fcidx = 1:length(flowcellfile); end

for j=1:length(fcidx), i = fcidx(j);
  fprintf('\nRunning SeqSeq for flowcell %s\n', flowcell{i});
  try
    SegSeq('-s',flowcell{i},'-t',flowcellfile{i});
  catch me
    fprintf('Problem in SegSeq\n');
    disp(me);
    disp(me.message);
  end
end

%cd([base_dir '/matfiles']);
%%
%
%foreach FLOWCELLFILE (*tumor_normal_aligned_paired_reads*mat)
%
%setenv FLOWCELL `echo $FLOWCELLFILE | sed 's/_tumor.*$//g' | sed 's/_medianfc//g'| sed 's/_nonmedfc//g'
%
%matlab -nodisplay -r "addpath ~gsaksena/CancerGenomeAnalysis/trunk/segseq; try SegSeq -s $FLOWCELL -t $FLOWCELLFILE
%...
%      -wig 10);catch, sin(1); end; exit;"
%
%end
