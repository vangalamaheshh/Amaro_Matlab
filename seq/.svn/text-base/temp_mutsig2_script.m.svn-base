function G = temp_mutsig2_script(P)

if ~exist('P','var'), P=[]; end

% required parameters
P = impose_default_value(P,'outdir','*required*');
P = impose_default_value(P,'mutfile','*required');

% data parameters
P = impose_default_value(P,'individual_set_id','mutsig');
P = impose_default_value(P,'banner','M:');
P = impose_default_value(P,'build','hg19');
P = impose_default_value(P,'build_dir',['/xchip/cga/annotation/db/ucsc/' P.build]);
P = impose_default_value(P,'terrfile','/cga/tcga-gsc/home/lawrence/mut/analysis/20110909_pancan/preprocess/pancan.terr_only.mat');
P = impose_default_value(P,'impute_full_coverage',true);
P = impose_default_value(P,'summed_cov_track','');
P = impose_default_value(P,'genes_to_analyze','');
P = impose_default_value(P,'ignore_flank_indels',false);

% computation parameters
P = impose_default_value(P,'interactive',false);
P = impose_default_value(P,'mutsig2_keyboard_every_report',false);  % (only applicable in interactive mode)
P = impose_default_value(P,'jobcount',100);    % (if not in interactive mode)
P = impose_default_value(P,'mutsig2_maxperm',1e8);  % (early iterations will start out lower)
P = impose_default_value(P,'mutsig2_single_iteration',P.mutsig2_maxperm<3e5);
P = impose_default_value(P,'bsub_mem',7);
P = impose_default_value(P,'bsub_queue_scatter','hour');
P = impose_default_value(P,'mutsig_version',[]);
P = impose_default_value(P,'automatically_add_zero_mutation_patients',false);
P = impose_default_value(P,'perform_mutsig2_analysis',true);
P = impose_default_value(P,'skip_directly_to_mutsig2_analysis',true);
P = impose_default_value(P,'mutsig2_maxsec',3600);  % (early iterations will start out lower)
P = impose_default_value(P,'allow_kill_long_jobs',true);
P = impose_default_value(P,'mutsig2_randseed',round(1e9*rand));
P = impose_default_value(P,'mutsig2_theta',1);
P = impose_default_value(P,'bsub_priority',99);
P = impose_default_value(P,'maxtries',2);

if P.interactive
  P = rmfield(P,'jobcount');
  ede(P.outdir);params_file = [P.outdir '/params.txt']; write_params_file(P,params_file);
  args={'-b',P.build,'-bd',P.build_dir,'-maf',P.mutfile,'-cov',P.terrfile,'-p',params_file};
  if isfield(P,'mutsig_version') && ~isempty(P.mutsig_version) && P.mutsig_version(1)=='/'
    libdir = P.mutsig_version
  else
    libdir = '/cga/tcga-gsc/home/lawrence/MutSigRun.203.5963';
    % libdir needed only for Java classes: does not need to be the most updated Matlab build
  end
  M = fh_MutSigRun(libdir,args{:})
  return
end

if isfield(P,'genelist') && ~isempty(P.genelist)
  error('for temp_mutsig2_script.m, please use P.genes_to_analyze, not P.genelist.');
end

if ~isempty(P.genes_to_analyze)
  if iscellstr(P.genes_to_analyze)
    genes_to_analyze = P.genes_to_analyze;
    P = rmfield(P,'genes_to_analyze');
  else
    error('invalid format for P.genes_to_analyze');
  end
end

if any(P.banner<32) error('P.banner contains unprintable characters'); end

iter = 1;
maxperm = 1000;
jobsec = 3600;  % each job should take <=1 hr

while(true)
  fprintf('Run %s  Iteration %d:   ',P.banner,iter);
  PP=P;
  PP.outdir = [P.outdir '/iter' num2str(iter)]; ede(PP.outdir);

  if ~exist('genes_to_analyze','var')
    fprintf('all genes');
    ngenes = 18000;
  else
    ngenes = length(genes_to_analyze);
    fprintf('%d genes',ngenes);
    PP.genelist = [PP.outdir '/genes_to_run.txt'];
    save_lines(genes_to_analyze,PP.genelist);
    PP.jobcount = min(PP.jobcount,ngenes);
  end

  genes_per_job = ceil(ngenes/PP.jobcount);
  if genes_per_job==1
    if ~strcmp(P.bsub_queue_scatter,'hour'), jobsec = jobsec * 2; end
    PP.min_before_kill = inf;
    maxperm = P.mutsig2_maxperm;  % if we're down to a single job, don't interrupt the permutations before max reached!
 else
    PP.mutsig2_maxsec = min(P.mutsig2_maxsec,ceil(20*jobsec/genes_per_job));  % allow 20x time factor for long genes
    if P.allow_kill_long_jobs
      PP.min_before_kill = round((jobsec*2)/60);   % jobs that run for more than 2x the max time are killed+restarted
    else
      PP.min_before_kill = inf;
    end
  end

  fprintf('   %d jobs   maxperm %d\n',PP.jobcount,maxperm);
  PP.banner = [P.banner num2str(iter) '.'];
  PP.mutsig2_maxperm = maxperm;
  PP.mutsig2_randseed = P.mutsig2_randseed + (iter-1);
  
  % scatter
  PP.wait_for_jobs = true;
  PP.skip_gather = true;
  result = bsub_mutsig_scatter_gather(PP);
  if result==-1
    fprintf('Note: some jobs failed.  Checking to see if all jobs produced usable output.\n');
    maxtryi=10;
    for tryi=1:maxtryi
      d=[]; d.file = direc([PP.outdir '/sc*/*.out']);
      if slength(d)<PP.jobcount
        fprintf('Found only %d/%d output files... waiting 1 min and trying again (try %d/%d)\n',...
                slength(d),PP.jobcount,tryi,maxtryi);
        pause(60);  % pause 1 minute and try again: maybe LSF is still writing output files
      else
        fprintf('Found all %d output files.\n',PP.jobcount);
        break
      end
    end
    d.bad = false(slength(d),1); d.good=d.bad;
    for i=1:slength(d)
      z = load_lines(d.file{i});
      if isempty(grep('REPORT',z))
        d.bad(i)=true;
      else
        d.good(i)=true;
      end
    end
    if sum(d.good)<PP.jobcount
      fprintf('%d jobs seem to have never produced any usable output.  Type CTRL-C to intervene.\n',PP.jobcount-sum(d.good));
      pause(600); % 10 min
      fprintf('Continuing.  Be aware some genes may have incomplete information.\n');
    else
      % maybe they just ran out of time
      fprintf('Looks OK... proceeding\n');
    end
  end

  if 0  % manually refresh old results
    num=5;
    G=cell(num,1);for i=1:num, tmp=load([P.outdir '/iter' num2str(i) '/results.mat']);G{i}=tmp.G;end
    G = concat_structs_keep_all_fields(G);
    G = rmfield(G,'line');
  end

  % gather
  newG = temp_mutsig2_gather(PP.outdir);
  if ~exist('G','var')
    G = newG;
  else
    % combine with previous iterations
    G = concat_structs_keep_all_fields({G,newG});
    [u ui uj] = unique(G.gene);
    keep = true(slength(G),1);
    for i=1:length(u), j=find(uj==i);
      G.nperm(j(1)) = sum(G.nperm(j));
      flds={'cons','clust','joint'}; for f=1:length(flds), fld=flds{f};
        G.(['k' fld])(j(1)) = sum(G.(['k' fld])(j));
      end
      keep(j(2:end))=false;
    end
    G = reorder_struct(G,keep);
    flds={'cons','clust','joint'}; for f=1:length(flds), fld=flds{f};
      [G.(['p' fld]) G.(['ci' fld])] = calc_pval_and_ci_ratio(G.(['k' fld]),G.nperm);
    end
  end

  % save analysis-so-far to the main directory (overwrite previous version)
  G = sort_struct(G,'pjoint');
  save_struct(G,[P.outdir '/results.txt']);
  save([P.outdir '/results.mat'],'G');

  % what genes need more permutations?
  G = sort_struct(G,'nmuts',-1); % to submit harder genes in earlier jobs
  gidx = find((G.cicons>P.mutsig2_theta | G.ciclust>P.mutsig2_theta | G.cijoint>P.mutsig2_theta) & (G.nperm<P.mutsig2_maxperm));
  if isempty(gidx), break; end
  genes_to_analyze = G.gene(gidx);
  maxperm = min(maxperm*10,P.mutsig2_maxperm);
  iter = iter + 1;

  % stopping condition is: all genes have ci<theta
  if P.mutsig2_single_iteration, break; end
end

fprintf('Done.\n');
