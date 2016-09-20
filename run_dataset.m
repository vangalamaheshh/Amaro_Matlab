function run_dataset(mat_fname,minvar,gamma,use_booster)

cd ~/projects/booster/datasets
load(mat_fname);
basename=mat_fname(1:(end-4));
c=cputime;

if use_booster
  bst='';
  [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr,topidx,K,N]=...
      marker_selection(D.dat,find(D.supdat(1,:)==min(D.supdat(1,:))),find(D.supdat(1,:)==max(D.supdat(1,:))),...
                       struct('method','ttest_minvar','minvar',minvar,'nparts_perm',1,...
                              'nparts_fix',1,'online',1,'booster',0.7,'gamma',gamma,'two_sided',1,...
                              'report',[basename ...
                      '_report.2side.10000.g95.txt']),10000);
else
  bst='no_boost.';
  [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr,topidx,K,N]=...
      marker_selection(D.dat,find(D.supdat(1,:)==min(D.supdat(1,:))),find(D.supdat(1,:)==max(D.supdat(1,:))),...
                       struct('method','ttest_minvar','minvar',minvar,'nparts_perm',1,...
                              'nparts_fix',1,'online',1,'report_th',0.7,'gamma',gamma,'two_sided',1,...
                              'report',[basename ...
                      '_report.2side.10000.g95.' bst 'txt']),10000);
end  
tot=cputime-c;
save([ basename '_res.2side.10000.g95.' bst 'mat'],'K','N','gamma','minvar','tot');
