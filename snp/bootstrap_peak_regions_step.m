function [regs,pvs]=bootstrap_peak_regions_step(CL,score_type,qv_thresh,ext)

[q,p,d,ads]=snp_score_permutations(CL,score_type,-1);
%  qv_thresh=0.25;

% plot_snp_score([ num2str(size(CL.dat,2)) '_bootstrp_' num2str(bi) ext],...
%                 CL,q,ads,qv_thresh,1,1); % 0
for k=1:2
  score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));
end
regs=generate_regs_by_peel_off(CL,d,q,score_type,score_thresh,501);
pvs=q;
