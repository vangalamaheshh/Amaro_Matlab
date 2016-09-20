function [D,gsupids]=add_pharma_scores(D,supids,score_type)

if ~exist('score_type','var')
  score_type=struct('method','ttest_minvar','minvar',0.4^2,'nparts_perm',1,'nparts_fix',1,'online',1,'booster', ...
                    0.7,'two_sided',1,'no_report','test_report.txt');
end

gsupids=[];
for i=1:length(supids)
  [idx,q,p,s,pi0,F]=get_top_markers(D,supids(i),...
                                      score_type,...
                                      -1,...
                                      struct('method','top','n',100,'thresh',0.25));
  
  [D,gsupids(end+1)]=add_D_sup(D,['SC-' deblank(D.supacc(supids(i),:))],...
                               [ 'Score for ' deblank(D.supdesc(supids(i),:))],F.s','rows');
  [D,tmp]=add_D_sup(D,['P-' deblank(D.supacc(supids(i),:))],[ 'P-value for ' deblank(D.supdesc(supids(i),:))],F.p','rows');
end

