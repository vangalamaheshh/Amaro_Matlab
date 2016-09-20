function RUN = temp_mutsig2_metascript(P,assume_bsub_empty_flag)
% same as temp_mutsig2_script
% but handles all runs at once: pancan, perttype, and leaveoneout
  
if ~exist('P','var'), P=[]; end

% required parameters
P = impose_default_value(P,'outdir','*required*');
P = impose_default_value(P,'mutfile','*required');

% data parameters
P = impose_default_value(P,'RUN',[]); % supersedes P.run_* and P.jobcount_* and P.bsub_mem_*
P = impose_default_value(P,'run_pancan',true);
P = impose_default_value(P,'run_perttype',true);
P = impose_default_value(P,'run_leaveoneout',true);
P = impose_default_value(P,'individual_set_id','mutsig');
P = impose_default_value(P,'build','hg19');
P = impose_default_value(P,'build_dir',['/xchip/cga1/annotation/db/ucsc/' P.build]);
P = impose_default_value(P,'terrfile',['/cga/tcga-gsc/home/lawrence/mut/analysis/20110909_pancan/preprocess/pancan.terr_only.' P.build '.mat']);
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
P = impose_default_value(P,'jobcount_pan',P.jobcount);
P = impose_default_value(P,'jobcount_ttype',P.jobcount);
P = impose_default_value(P,'bsub_mem_pan',P.bsub_mem);
P = impose_default_value(P,'bsub_mem_ttype',P.bsub_mem);
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
P = impose_default_value(P,'add_genes_missing_from_iter_results',true); % in case of aborted jobs.  but still getting stuck on nmut<2 genes

ede(P.outdir);

if strcmp(P.build,'hg18') && grepm('/cga/tcga-gsc/home/lawrence/mut/analysis/20121031_pancan/cbb/exponential',{P.summed_cov_track})
  error('Looks like you''re using an hg19 summed_cov_track with hg18 mutations');
end

if P.interactive
  fprintf('************************************\n');
  fprintf('*                                  *\n');
  fprintf('*   STARTING IN INTERACTIVE MODE   *\n');
  fprintf('*                                  *\n');
  fprintf('************************************\n');
  P = rmfield(P,'jobcount');
  params_file = [P.outdir '/params.txt']; write_params_file(P,params_file);
  args={'-b',P.build,'-bd',P.build_dir,'-maf',P.mutfile,'-cov',P.terrfile,'-p',params_file};
  if isfield(P,'mutsig_version') && ~isempty(P.mutsig_version) && P.mutsig_version(1)=='/'
    libdir = P.mutsig_version
  else
    libdir = '/cga/tcga-gsc/home/lawrence/MutSigRun.6.5971';
    % libdir needed only for Java classes: does not need to be the most updated Matlab build
  end
  M = fh_MutSigRun(libdir,args{:})
  return
end

run_results_file = [P.outdir '/results.mat'];

if exist(run_results_file,'file')
  % pick up from where we previously left off
  fprintf('Resuming...\n');
  load(run_results_file,'RUN');
end

if ~exist('RUN','var')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PRE-ANALYZE MAF FILE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Loading maf file...\n');
  M=[]; M.mut = load_struct(P.mutfile);
  if P.run_perttype || P.run_leaveoneout || ~isempty(P.RUN)
    if ~isfield(M.mut,'ttype'), error('need "ttype" field in MAF in order to run per_ttype or leaveoneout analysis'); end
    [M.ttype.name tmp M.mut.ttype_idx] = unique(M.mut.ttype);
    if ~isempty(P.RUN)
      save_struct(P.RUN,[P.outdir '/runs.txt']);
    else
      save_struct(M.ttype,[P.outdir '/ttypes.txt']);
    end
  end
  M.mut = demand_field_with_replace(M.mut,'gene','Hugo_Symbol');
  M.mut = reorder_struct_exclude(M.mut,strcmp('?',M.mut.gene)|strcmp('',M.mut.gene)|strcmpi('unknown',M.mut.gene));
  M.mut = reorder_struct_exclude(M.mut,grepmi('synon|silent|UTR|flank|intron|IGR',M.mut.type));  % not used in MutSig2
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
  if exist('genes_to_analyze','var')
    M.mut = reorder_struct(M.mut,ismember(M.mut.gene,genes_to_analyze));
  end
  % filter against master genelist from terr file
  tmp = load(P.terrfile,'C1');
  M.mut = reorder_struct(M.mut,ismember(M.mut.gene,tmp.C1.gene.name));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PLAN THE RUNS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Planning the runs, writing MAFs for each run...\n');

  if isempty(P.RUN)
    RUN = {}; ri=0;
    if P.run_pancan
      ri=ri+1;
      RUN.name{ri,1} = 'pancan';
      RUN.ttype{ri,1} = 'pancan';
      RUN.banner{ri,1} = 'PAN';
      RUN.outdir{ri,1} = [P.outdir '/' RUN.name{ri}]; ede(RUN.outdir{ri});
      RUN.mutfile{ri,1} = [RUN.outdir{ri} '/' RUN.name{ri} '.maf'];
      if ~exist(RUN.mutfile{ri},'file') save_struct(M.mut,RUN.mutfile{ri}); end
      RUN.genes{ri,1} = unique(M.mut.gene);
      RUN.jobcount(ri,1) = P.jobcount_pan;
      RUN.bsub_mem(ri,1) = P.bsub_mem_pan; 
    end
    if P.run_perttype
      for tti=1:slength(M.ttype)
        ri=ri+1;
        RUN.name{ri,1} = ['ttype' num2str(tti)];
        RUN.ttype{ri,1} = M.ttype.name{tti};
        RUN.banner{ri,1} = ['TT' num2str(tti)];
        RUN.outdir{ri,1} = [P.outdir '/' RUN.name{ri}]; ede(RUN.outdir{ri});
        RUN.mutfile{ri,1} = [RUN.outdir{ri} '/' RUN.name{ri} '.maf'];
        midx = find(M.mut.ttype_idx==tti);
        if ~exist(RUN.mutfile{ri},'file') save_struct(reorder_struct(M.mut,midx),RUN.mutfile{ri}); end
        RUN.genes{ri,1} = unique(M.mut.gene(midx));
        RUN.jobcount(ri,1) = P.jobcount_ttype;
        RUN.bsub_mem(ri,1) = P.bsub_mem_ttype;
    end,end
    if P.run_leaveoneout
      for tti=1:slength(M.ttype)
        ri=ri+1;
        RUN.name{ri,1} = ['leaveout' num2str(tti)];
        RUN.ttype{ri,1} = ['leaveout_' M.ttype.name{tti}];
        RUN.banner{ri,1} = ['LO' num2str(tti)];
        RUN.outdir{ri,1} = [P.outdir '/' RUN.name{ri}]; ede(RUN.outdir{ri});
        RUN.mutfile{ri,1} = [RUN.outdir{ri} '/' RUN.name{ri} '.maf'];
        midx = find(M.mut.ttype_idx~=tti);
        if ~exist(RUN.mutfile{ri},'file') save_struct(reorder_struct(M.mut,midx),RUN.mutfile{ri}); end
        RUN.genes{ri,1} = unique(M.mut.gene(midx));
        RUN.jobcount(ri,1) = P.jobcount_pan;
        RUN.bsub_mem(ri,1) = P.bsub_mem_pan;
    end,end
    if slength(RUN)==0, error('Nothing to do!'); end

  else  % P.RUN was provided
    fprintf('(using P.RUN)\n');
    RUN = P.RUN;
    demand_fields(RUN,{'name','outdir','jobcount','bsub_mem','ttype_grep','ttype_grepv'});
    if any(grepm('^/',RUN.outdir)), error('P.RUN.outdir should be relative path'); end
    RUN.ttype = RUN.name;
    RUN.banner = str2cell(sprintf('R%02d\n',1:slength(RUN)));
    for ri=1:slength(RUN)
      RUN.outdir{ri,1} = [P.outdir '/' RUN.outdir{ri}]; ede(RUN.outdir{ri});
      RUN.mutfile{ri,1} = [RUN.outdir{ri} '/' RUN.name{ri} '.maf'];
      midx = (1:slength(M.mut))';
      if ~isempty(RUN.ttype_grep{ri}), midx = midx(grepm(RUN.ttype_grep{ri},M.mut.ttype(midx))); end
      if ~isempty(RUN.ttype_grepv{ri}), midx = midx(~grepm(RUN.ttype_grepv{ri},M.mut.ttype(midx))); end
      if ~exist(RUN.mutfile{ri},'file') save_struct(reorder_struct(M.mut,midx),RUN.mutfile{ri}); end
      RUN.genes{ri,1} = unique(M.mut.gene(midx));
    end
    if slength(RUN)==0, error('Nothing to do!'); end
  end

  % remove jobs with no genes to run
  RUN = reorder_struct_exclude(RUN,cellfun('length',RUN.genes)==0);

  RUN.iter = ones(slength(RUN),1);
  RUN.maxperm = 1000*ones(slength(RUN),1);
  RUN.jobsec = 3600*ones(slength(RUN),1);    % initially each job should take <=1 hr; note, this parameter is not actually being used!
  RUN.genes_remaining = RUN.genes;
  RUN.G = cell(slength(RUN),1);
  RUN.finished = false(slength(RUN),1);

  save(run_results_file,'RUN');
end

P = rmfield_if_exist(P,'RUN'); % (causes problems when saving to params.txt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCATTER/GATHER THE RUNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RUN.first_time_examining_run_during_this_session = true(slength(RUN),1);

do_this_run_first = nan;
all_finished = false;
while(~all_finished)

  % sort the RUNS so that we process them in order of 
  %            earliest iteration to latest iteration
  % (but make an exception for run designated by do_this_run_first)
  z = keep_fields(RUN,{'iter','finished'});
  z.ri = as_column(1:slength(z));
  z = sort_struct(z,'iter');
  z.is_first = 1*(z.ri==do_this_run_first);
  z = sort_struct(z,'is_first',-1);
  run_order = z.ri;

  % iterate over the runs
  % 1. see if the current iteration already has results on disk
  %    -- if so, re-load those results and go to the next iteration
  % 2. if not, if the current iteration is still running
  %    -- if so, see if all jobs have ended
  %    -- if not, see if some have gone over time
  %    -- if not, we're done.
  % 3. if all the jobs have ended, then read them in and write the results,
  %    -- update RUN
  %    -- go on to the next iteration
  % 4. if the current iteration has no jobs running, then submit
  
  all_finished = true;
  for ridx=1:slength(RUN), ri=run_order(ridx);
    if RUN.finished(ri), continue; end
    all_finished = false;

    iter = RUN.iter(ri);
    fprintf('---------%d/%d (%s) %s  Iter %d',ri,slength(RUN),RUN.banner{ri},RUN.ttype{ri},iter);
    iterdir = [RUN.outdir{ri} '/iter' num2str(iter)]; ede(iterdir);
    iter_results_file = [iterdir '/results.mat'];
    if exist(iter_results_file,'file')
      fprintf('  retrieving saved results\n');
      newG = load(iter_results_file,'G'); newG=newG.G;
      iteration_finished = true;
    else      
      % results file doesn't exist.  so this means we need to:
      % call the main bsub script, which will 
      % see if jobs need to be submitted or resubmitted (if so, submits them and returns result=2)
      % or if we're just waiting (result=2) 
      % or if it's time to gather (result=0,all succeeded or result=-1,some failures)
      PP=P;
      PP.outdir = iterdir;
      PP.mutfile = RUN.mutfile{ri};
      genes_to_analyze = RUN.genes_remaining{ri};
      ngenes = length(genes_to_analyze);
      PP.genelist = [PP.outdir '/genes_to_run.txt'];
      save_lines(genes_to_analyze,PP.genelist);
      PP.jobcount = min(RUN.jobcount(ri),ngenes);
      genes_per_job = ceil(ngenes/PP.jobcount);
      if genes_per_job==1
        if ~strcmp(P.bsub_queue_scatter,'hour'), RUN.jobsec(ri) = RUN.jobsec(ri) * 2; end
        PP.min_before_kill = inf;
        RUN.maxperm(ri) = P.mutsig2_maxperm;  % if we're down to a single job, don't interrupt the permutations before max reached!
      else
        PP.mutsig2_maxsec = min(P.mutsig2_maxsec,ceil(20*RUN.jobsec(ri)/genes_per_job));  % allow 20x time factor for long genes
        if P.allow_kill_long_jobs
          PP.min_before_kill = round((RUN.jobsec(ri)*2)/60);   % jobs that run for more than 2x the max time are killed+restarted
        else
          PP.min_before_kill = inf;
        end
      end
      PP.banner = [RUN.banner{ri} ':'  num2str(iter) '.'];
      PP.bsub_mem = RUN.bsub_mem(ri);
      PP.mutsig2_maxperm = RUN.maxperm(ri);
      PP.mutsig2_randseed = P.mutsig2_randseed + (iter-1);
      PP.wait_for_jobs = false;
      PP.skip_gather = true;
      if exist('assume_bsub_empty_flag','var') && assume_bsub_empty_flag && RUN.first_time_examining_run_during_this_session(ri)
        PP.assume_no_jobs_running_or_pending = true;
      end
      % CALL submit script
      fprintf('  %d genes  %d jobs @ maxperm %d  ',ngenes,PP.jobcount,PP.mutsig2_maxperm);
      if PP.jobcount>0
        result = bsub_mutsig_scatter_gather(PP);
      else
        fprintf('WARNING:  jobcount==0\n');
        result = 2;
      end
      RUN.first_time_examining_run_during_this_session(ri) = false;
      if result==2
        % this means all jobs have been submitted and we're now waiting
        iteration_finished = false;
      else
        % else it's probably time to gather: first check if we need to warn about failed jobs
        if result==-1   % some jobs failed
          fprintf('Checking to see if all jobs produced usable output.\n');
          maxtryi=10;
          for tryi=1:maxtryi
            d=[]; d.file = direc([PP.outdir '/sc*/*.out']);
            if slength(d)<PP.jobcount
              fprintf('Found only %d/%d output files... waiting 1 min and trying again (try %d/%d)\n',...
                      slength(d),PP.jobcount,tryi,maxtryi);
              pause(60);  % pause 1 minute and try again: maybe LSF is still writing output files
            else
              fprintf('Found all %d output files.\n',PP.jobcount); break
            end
          end
          d.bad = false(slength(d),1); d.good=d.bad;
          for i=1:slength(d)
            z = load_lines(d.file{i});
            if isempty(grep('REPORT',z)), d.bad(i)=true; else d.good(i)=true; end
          end
          if sum(d.good)<PP.jobcount
            fprintf('%d jobs seem to have never produced any usable output.\n',PP.jobcount-sum(d.good));
          else
            fprintf('Looks OK... proceeding\n');   % maybe they just ran out of time
          end
        end
        % GATHER 
        newG = temp_mutsig2_gather(PP.outdir);
        iteration_finished = true;
      end
    end

    if iteration_finished
      if iter==1
        % this is the first iteration: use newG directly
        G = newG;
        if ~isfield(G,'pjoint'), error('PROBLEM! ALL JOBS IN ITER1 FAILED!'); end
      else
        % else, combine with previous iterations
        G = RUN.G{ri};
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
      G = sort_struct(G,{'pjoint','nmuts'},[1 -1]);
      RUN.G{ri} = G;
      % save analysis-so-far to the main directory (overwrite previous version)
      G.q = calc_fdr_value(G.pjoint);
      save_struct(G,[RUN.outdir{ri} '/results.txt']);
      save([RUN.outdir{ri} '/results.mat'],'G');
      % what genes need more permutations?
      gidx = find((G.cicons>P.mutsig2_theta | G.ciclust>P.mutsig2_theta | G.cijoint>P.mutsig2_theta) & (G.nperm<P.mutsig2_maxperm));
      genes_needing_more = G.gene(gidx);
      missing_genes = setdiff(RUN.genes{ri},G.gene);
      % stopping condition is: all genes have ci<theta
      if isempty(genes_needing_more) || (P.mutsig2_single_iteration && isempty(missing_genes))
        RUN.finished(ri) = true;
      else
        % before starting next run: any genes missing entirely? add them (unless otherwise specified)
        if ~isempty(missing_genes)
          if P.add_genes_missing_from_iter_results
            fprintf('Note: adding %d missing genes to the next iteration\n',length(missing_genes));
          else
            fprintf('%d genes are missing, but not adding these to the next iteration.  (maybe they are nmut<2)\n',length(missing_genes));
            missing_genes=[];
          end
        end
        % update RUN info and increment iteration conuter
        RUN.genes_remaining{ri} = [genes_needing_more;missing_genes];
        if length(RUN.genes_remaining{ri})<8000
          RUN.maxperm(ri) = min(RUN.maxperm(ri)*10,P.mutsig2_maxperm);
        end
        RUN.iter(ri) = RUN.iter(ri) + 1;
      end
      % save updated run results
      save([run_results_file '.partial'],'RUN');
      system(['mv ' run_results_file '.partial ' run_results_file]);

      do_this_run_first = ri;
      break    % restart scan and submit this run's new jobs
    end % if iteration_finished
    
  end   % examine next run
end     % while(~all_finished)

fprintf('Done.\n');


