function R = mini_driver(M, P)
%The main function for the clustering program.

if ~iscell(M), M={M}; end

if ischar(P) || isnumeric(P)
  tmp=P;
  P=[];
  P.gene=tmp;
end

if ~exist('P','var'), P=[]; end

P = impose_default_value(P, 'gene', (1:slength(M{1}.gene)));
P = impose_default_value(P, 'setnames', repmat({'iset'},length(M),1));
P = impose_default_value(P, 'mutsig2_randseed',1234);
P = impose_default_value(P, 'mutsig2_restrict_to_nonsilent_coding',true);
P = impose_default_value(P, 'mutsig2_restrict_to_one_mutation_per_patient',false); % should be set to true!!!!!
P = impose_default_value(P, 'conservation_hg19_fwb','/cga/tcga-gsc/home/lawrence/db/hg19/conservation46/all.fwb');
P = impose_default_value(P, 'conservation_hg18_fwb','/cga/tcga-gsc/home/lawrence/db/hg18/conservation44/all.fwb');
P = impose_default_value(P, 'context_and_sense_hg18_fwb','/cga/tcga-gsc/home/lawrence/db/hg18/c65e29/all.fwb');
P = impose_default_value(P, 'context_and_sense_hg19_fwb','/cga/tcga-gsc/home/lawrence/db/hg19/c65e29/all.fwb');
P = impose_default_value(P, 'context_and_sense_hg18_categs','/cga/tcga-gsc/home/lawrence/db/hg18/c65e29/categs.txt');
P = impose_default_value(P, 'context_and_sense_hg19_categs','/cga/tcga-gsc/home/lawrence/db/hg19/c65e29/categs.txt');
P = impose_default_value(P, 'mutsig2_use_sample_specific_coverage',false);
P = impose_default_value(P, 'mutsig2_use_power_calculation', false);
P = impose_default_value(P, 'mutsig2_power_calculation_method', 'lod_score');
P = impose_default_value(P, 'impute_full_coverage',false);
P = impose_default_value(P, 'radius_to_impute_coverage_around_mutations',10);
P = impose_default_value(P, 'min_num_samples_to_median_for_unknown_coverage',100);
P = impose_default_value(P, 'mutsig2_remove_nan_allelic_frac', false);

%% Make sure power calculation and sample specific coverage are not both on at the same time 

if P.mutsig2_use_sample_specific_coverage && P.mutsig2_use_power_calculation
  error('"Power calculation" and "sample specific coverage" are two different methods, cannot use both');
end 

%%% Make sure sample specific coverage requirements are met 
if P.mutsig2_use_sample_specific_coverage
  allok = true;
  if ~isfield(P,'mutsig2_coverage_fwi')
    fprintf('Need P.mutsig2_coverage_fwi to use sample-specific coverage\n'); allok = false;
  else
    if ~exist(P.mutsig2_coverage_fwi,'file')
      fprintf('Not found: P.mutsig2_coverage_fwi = %s\n',P.mutsig2_coverage_fwi); allok = false;
    end
    for si=1:length(M)
      if ~isfield(M{si},'pat') && isfield(M{si},'patient'), M{si}=rename_field(M{si},'patient','pat'); end
      if ~isfield(M{si},'pat')
        fprintf('Need M.pat to use sample-specific coverage\n'); allok=false; break;
      end
        if ~isfield(M{si}.mut,'patient')
          fprintf('Need M.mut.patient to use sample-specific coverage\n'); allok=false; break;
        end
        if ~isfield(M{si}.pat,'name')
          fprintf('Need M.pat.name to use sample-specific coverage\n'); allok=false; break;
        end          
        M{si}.mut.pat_idx = listmap(M{si}.mut.patient,M{si}.pat.name);
        if all(isnan(M{si}.mut.pat_idx))
          fprintf('M.mut.pat_idx is bad: can''t use sample-specific coverage\n'); allok=false; break;
        end
        if ~isfield(M{si}.pat,'fwb')
          fprintf('Need M.pat.fwb to use sample-specific coverage\n'); allok=false; break;
        end
        fprintf('Verifying presence of sample-specific coverage files:\n');
        M{si}.pat.fwb_ok = demand_files(M{si}.pat.fwb);
        nok = sum(M{si}.pat.fwb_ok);
        if nok==0
          fprintf('No M.pat.fwb exists!  Cannot use sample-specific coverage\n'); allok=false; break;
  end,end,end
  if ~allok
    fprintf('***WARNING***   Setting P.mutsig2_use_sample_specific_coverage = false\n');
    P.mutsig2_use_sample_specific_coverage = false;
  end
end


%%% Make sure power calculation requirements are met 
%% both sample-specific coverage and power calculation will use two different 


if P.mutsig2_use_power_calculation
  allok = true;
  if ~isfield(P,'mutsig2_power_coverage_fwi')
    fprintf('Need P.mutsig2_power_coverage_fwi to use power calculation method\n'); allok = false;
  else
    if ~exist(P.mutsig2_power_coverage_fwi,'file')
      fprintf('Not found: P.mutsig2_power_coverage_fwi = %s\n',P.mutsig2_coverage_fwi); allok = false;
    end
    for si=1:length(M)
      if ~isfield(M{si},'pat') && isfield(M{si},'patient'), M{si}=rename_field(M{si},'patient','pat'); end
      if ~isfield(M{si},'pat')
        fprintf('Need M.pat to use power calculation\n'); allok=false; break;
      end
        if ~isfield(M{si}.mut,'patient')
          fprintf('Need M.mut.patient to use power calculation\n'); allok=false; break;
        end
        if ~isfield(M{si}.pat,'name')
          fprintf('Need M.pat.name to use power calculation'); allok=false; break;
        end
        M{si}.mut.pat_idx = listmap(M{si}.mut.patient,M{si}.pat.name);
        if all(isnan(M{si}.mut.pat_idx))
          fprintf('M.mut.pat_idx is bad: can''t use power calculation\n'); allok=false; break;
        end
        if ~isfield(M{si}.pat,'power_files')
          fprintf('Need M.pat.power_files to use power calculation\n'); allok=false; break;
        end
        fprintf('Verifying presence of sample-specific coverage files:\n');
        M{si}.pat.power_ok = demand_files(M{si}.pat.power_files);
        nok = sum(M{si}.pat.power_ok);
        if nok==0
          fprintf('No M.pat.power_files exists!  Cannot use sample-specific coverage\n'); allok=false; break;
  end,end,end
  if ~allok
    fprintf('***WARNING***   Setting P.mutsig2_use_power_calculation = false\n');
    P.mutsig2_use_power_calculation = false;
  end
end




if P.mutsig2_use_sample_specific_coverage && P.impute_full_coverage
  fprintf('Not using full coverage, because using sample-specific coverage\n');
  P.impute_full_coverage = false;
end

if P.mutsig2_use_power_calculation && P.impute_full_coverage
  fprintf('Not using full coverage, because using power calculation\n');
  P.impute_full_coverage = false;
end



if ~P.mutsig2_use_sample_specific_coverage   % use summed coverage profile
  if ~P.impute_full_coverage
    allcovpresent = true;
    for si=1:length(M)
      if ~isfield(M{si},'file') || ~isfield(M{si}.file,'summed_cov_track') || isempty(M{si}.file.summed_cov_track)
        allcovpresent = false;
      end
    end
    if ~allcovpresent
      fprintf('Missing summed_cov_track(s): will force P.impute_full_coverage\n');
      P.impute_full_coverage = true;
    end
  end
  P = impose_default_value(P, 'mutsig2_max_coverage_bins',20);   % (see below for explanation)
  if P.impute_full_coverage, P.mutsig2_max_coverage_bins = 1; end
end

if ~isnumeric(P.gene)
    tmp = listmap(P.gene,M{1}.gene.name);
    if all(isnan(tmp)), error('invalid P.gene'); end
    P.gene = tmp(~isnan(tmp));
end

nsets = length(M);
builds = cell(nsets, 1);

% INITIALIZE RANDOM NUMBER GENERATOR

if P.mutsig2_randseed<1
  fprintf('P.mutsig2_randseed<1:  resetting to a random integer value\n');
  P.mutsig2_randseed = round(1e9 * rand());
end
fprintf('MutSig2: initializing with randseed=%d\n',P.mutsig2_randseed);
rand('twister',P.mutsig2_randseed);

% PREPROCESS EACH DATASET

for si=1:nsets
  builds{si} = M{si}.build;
  M{si}.mut.chr = convert_chr(M{si}.mut.chr);
  M{si}.mut = make_numeric(M{si}.mut,{'start','end','categ', 'i_tumor_f'});

  % NOTE: copied from call_mutsig_cluster.m
  if P.mutsig2_restrict_to_nonsilent_coding
    fprintf('Restricting to nonsilent coding mutations\n');
    if isfield(M{si},'use_nonsilent')
      % keep only use_nonsilent mutations
      M{si}.mut = reorder_struct(M{si}.mut,M{si}.use_nonsilent);
      M{si} = rmfield(M{si},'use_nonsilent');
    elseif isfield(M{si}.mut,'is_silent') && isfield(M{si}.mut,'is_coding')
      M{si}.mut = reorder_struct(M{si}.mut,M{si}.mut.is_coding & ~M{si}.mut.is_silent);
    elseif isfield(M{si}.mut,'type')
      idx = grepi('missense|nonsense|splice|non.?stop|in.?frame|frame.?shift|de.?novo',M{si}.mut.type);
      M{si}.mut = reorder_struct(M{si}.mut,idx);
    else
      error('Don''t know how to select nonsilent coding mutations!');
    end
  end

end

hg19_idx = find(strcmp(builds, 'hg19'), 1,'first');
hg18_idx = find(strcmp(builds, 'hg18'), 1,'first');

%% Figure out what builds are present for the datasets used,
% and prepare a bi indexing variable to loop over builds for each gene
% for the coordinate processing bellow

if isempty(hg18_idx) && isempty(hg19_idx)
    error('No build information available!  Please set M.build / M{i}.build');
elseif isempty(hg18_idx)  % hg19 only
    bi = 1;
elseif isempty(hg19_idx)  % hg18 only
    bi = 2;
else
    bi = [1, 2];          % both hg18 and hg19
end

numgenes = length(P.gene);

R = [];
R.gene_name =  M{1}.gene.name(P.gene);
R.nperm = nan(numgenes,1);
R.p_clust = nan(numgenes,1);
R.p_cons = nan(numgenes,1);
R.p_joint = nan(numgenes,1);
R.eff_clust = nan(numgenes,1);
R.eff_cons = nan(numgenes,1);

if ~exist('first_idx','var'), first_idx = 1; end
if ~exist('last_idx','var'), last_idx = inf; end

ridx = 1;
for i = P.gene

    if i<first_idx || i>last_idx, continue; end
    
    gname = M{1}.gene.name{i};
    fprintf('%s\n',gname);
    
    % GENE characteristics
    starts = cell(2,1);
    ends = cell(2, 1);
    exons_matrix = cell(2, 1);
    genomic_array = cell(2,1);

    % CONSERVATION TRACKS
    file_conservation_hg19 = org.broadinstitute.cga.tools.seq.FixedWidthBinary(P.conservation_hg19_fwb);
    file_conservation_hg18 = org.broadinstitute.cga.tools.seq.FixedWidthBinary(P.conservation_hg18_fwb);
    
    % CONTEXT_AND_SENSE TRACKS
    file_context_and_sense_hg19 = org.broadinstitute.cga.tools.seq.FixedWidthBinary(P.context_and_sense_hg19_fwb);
    file_context_and_sense_hg18 = org.broadinstitute.cga.tools.seq.FixedWidthBinary(P.context_and_sense_hg18_fwb);
    context_and_sense_hg19_categs = load_struct(P.context_and_sense_hg19_categs);
    context_and_sense_hg18_categs = load_struct(P.context_and_sense_hg18_categs);
    
    error_found = 0;
    %% For each build, figure out the exons, the genomic coordinates as well as the polyphen and categories tracks %%
    fprintf('Processing conservation and coverage...\n');
    for bi = bi
      if bi == 1
        curr_idx = hg19_idx;
      else
        curr_idx = hg18_idx;
      end

      % find the mutations in this gene
      mut_chrs = cell(nsets,1);
      for si=1:nsets
        if isfield(M{si}.mut,'gene_idx')
          midx = find(M{si}.mut.gene_idx == i);
        else
          midx = find(M{si}.mut.gene == i);
        end
        mut_chrs{si} = convert_chr(M{si}.mut.chr(midx));
      end
      mut_chrs = cat(1,mut_chrs{:});
      mut_chr = unique(mut_chrs);

      % find out the chromosome of the mutations
      if length(mut_chr)>1
        fprintf('WARNING: gene %s chromosome assignment ambiguous in mutation list.\n',gname);
      end

      % find out the target list for this gene
      tidx = find(M{curr_idx(1)}.cov.targ.gidx == i);
      targ_chrs = M{curr_idx(1)}.cov.targ.chr(tidx);
      targ_chr = unique(targ_chrs);
      if length(targ_chr)>1
        fprintf('WARNING: gene %s chromosome assignment ambiguous in target list.\n',gname);
      end

      % handle edge case: pick a way to resolve chromosome ambiguity
      if length(mut_chr)>1 || length(targ_chr)>1
        isect_chr = intersect(mut_chr,targ_chr);
        if length(isect_chr)==0
          fprintf('No possible resolution of chromosome ambiguities for gene %s: skipping this gene\n',gname);
          fprintf('[REPORT]  %-11s   SKIP: ambiguous chromosomes\n',gname);
          continue
        elseif length(isect_chr)>1
          isect_ct = zeros(length(isect_chr),1);
          for j=1:length(isect_chr)
            isect_ct(j) = sum(mut_chrs==isect_chr(j));
          end
          [tmp idx] = max(isect_ct);
          isect_chr = isect_chr(idx);
        end
        fprintf('Choosing chr%d for analysis of gene %s\n',isect_chr,gname);
        idx = find(targ_chrs==isect_chr);
        tidx = tidx(idx);
        targ_chrs = targ_chrs(idx);
        targ_chr = isect_chr;
        idx = find(mut_chrs==isect_chr);
        mut_chrs = mut_chrs(idx);
        mut_chr = isect_chr;
      end          

      if isempty(tidx), continue; end
      
      % finish getting info for this gene
      starts{bi} = M{curr_idx(1)}.cov.targ.start(tidx);
      ends{bi} = M{curr_idx(1)}.cov.targ.end(tidx);
      curr_starts = starts{bi};
      curr_ends = ends{bi};
      genelength = 0;
      for j = 1:length(curr_starts)
        genelength = genelength + (curr_ends(j) - curr_starts(j)) + 1;
      end
      exons_matrix{bi} = [curr_starts'; curr_ends']';
      genomic_array{bi} = get_GenomicArray(exons_matrix{bi}, genelength);
      chr = mut_chr;
      if isempty(chr), chr=nan; end  % (allows graceful exit later)

      if bi == 1
        conservation_track{bi} = populate_conservation(genomic_array{bi}, file_conservation_hg19, chr);
        context_and_sense_track{bi} = populate_coverage(genomic_array{bi}, file_context_and_sense_hg19, chr);
        context_and_sense_categs{bi} = context_and_sense_hg19_categs;
      else
        conservation_track{bi} = populate_conservation(genomic_array{bi}, file_conservation_hg18, chr);
        context_and_sense_track{bi} = populate_coverage(genomic_array{bi}, file_context_and_sense_hg18, chr);
        context_and_sense_categs{bi} = context_and_sense_hg18_categs;
      end
    end  % next build
    
    % close files
    file_conservation_hg19.close();
    file_conservation_hg18.close();
    file_context_and_sense_hg19.close();
    file_context_and_sense_hg18.close();
    
    if error_found
       fprintf('Error with target list retrieved for gene: %s\n', gname);
       continue
    end
    
    %% Check the target list for errors (for sh-RNAs and such, mutliple target lists may come up as one)
    %error_found = 0;
    for k = 1:length(exons_matrix)
      for j = bi
        if ~isempty(starts{bi}) && (starts{bi}(1) > starts{bi}(length(starts{bi})))
          error_found = 1;
          break
        end
      end
    end
    if error_found
      fprintf('Error with target list retrieved for gene: %s\n', gname);
      continue
    end
    
    % MUTATIONS in each tumor type
    fprintf('Processing mutations... \n');
    mutpos = cell(nsets,1);
    patients = cell(nsets, 1);
    newbase = cell(nsets,1);
    num_mutations = zeros(nsets,1);
    mut_categories = cell(nsets,1);       % from nucleotide context
    mut_senses = cell(nsets,1);           % 0=indel  1=synon  2=missense  3=nonsense/readthrough/splice
    throw_flavor = cell(nsets,1);         % combines information of category + sense
    tot_num_mutations = 0;
    
    coverage_track = cell(nsets,1);
    coverage_track_factor = ones(nsets,1);
    maxCov = zeros(nsets,1);
    error_found = 0;
    skip_this_gene = false;
    for si = 1:nsets

      if strcmp(builds{si}, 'hg19')
        curr_exons_matrix = exons_matrix{1};
        curr_genomic_array = genomic_array{1};
      else
        curr_exons_matrix = exons_matrix{2};
        curr_genomic_array = genomic_array{2};
      end

      % find mutations in this gene
      if isfield(M{si}.mut,'gene_idx')
        midx = find(M{si}.mut.gene_idx == i);
      else
        midx = find(M{si}.mut.gene == i);
      end

      fprintf('Mutation breakdown by type:\n');
      count(M{si}.mut.type(midx));

      % see if any mutations in this territory are missing (for example, because annotated to another gene)
      % if so, add them to the list of mutations to permute
      midx2 = find(M{si}.mut.chr==chr & M{si}.mut.start<=curr_exons_matrix(end,end) & M{si}.mut.end>=curr_exons_matrix(1,1));
      midx2 = setdiff(midx2,midx);
      if ~isempty(midx2)
        fprintf('Including %d additional mutations that were not originally annotated to this gene.\n',length(midx2));
        midx = union(midx,midx2);
        fprintf('New mutation list:\n');
        xcount(M{si}.mut.type(midx),M{si}.mut.gene(midx));
      end

      % remove mutations that fail to map to an exon
      mutpos{si} = mutation_converter(curr_exons_matrix, M{si}.mut.start(midx));
      not_mapped = isnan(mutpos{si});
      if any(not_mapped)
        fprintf('Skipping %d mutations outside target territory of this gene.\n',sum(not_mapped));
        midx = midx(~not_mapped);
        mutpos{si} = mutpos{si}(~not_mapped);
        fprintf('Mutations remaining:\n');
        count(M{si}.mut.type(midx));
      end

      % if flag is up, remove mutations with allelic fraction = NaN 
      if P.mutsig2_remove_nan_allelic_frac 
        nan_AF = isnan(M{si}.mut.i_tumor_f(midx)); 
        midx = midx(~nan_AF);
        fprintf('Removing %d mutations with NaN allelic fraction.\n', sum(nan_AF));
        count(M{si}.mut.type(midx))
      end 

      % get rid of mutations with categ==0
      badcateg = M{si}.mut.categ(midx)==0 | isnan(M{si}.mut.categ(midx)) | M{si}.mut.categ(midx)>slength(M{si}.categ);
      if any(badcateg)
        fprintf('Skipping %d mutations outside the category set.\n',sum(badcateg));
        midx = midx(~badcateg);
        mutpos{si} = mutpos{si}(~badcateg);
        fprintf('Mutations remaining:\n');
        count(M{si}.mut.type(midx));
      end

      patients{si} = M{si}.mut.patient(midx);

      % keep only one (randomly chosen) mutation per patient %%%%%%%%%%%%%%%%%%
      if P.mutsig2_restrict_to_one_mutation_per_patient
        [upat ui upatj] = unique(patients{si});
        % randomize which one is kept (to avoid spurious clustering at the beginning of the gene)
        keep = ones(length(upat), 1);
        for qq=1:length(upat)
          idx = find(upatj==qq);
          if length(idx)>1
            tmp = randperm(length(idx));
            keep(qq) = tmp(1);
          end
        end
        fprintf('Keeping %d/%d mutations after restricting to one per patient\n',length(keep),length(upatj));
        muts_to_keep = nan(length(upat), 1);
        for k = 1:length(upat)
          if ~iscell(patients{si})
            pat_muts_idx = find(patients{si} == upat(k));
          else
            pat_muts_idx = find(strcmp(patients{si}, upat{k}));
          end
          pat_muts = mutpos{si}(pat_muts_idx); 
          muts_to_keep(k) = pat_muts(keep(k));
        end
        midx = midx(keep);
        mutpos{si} = muts_to_keep;
      end

      num_mutations(si) = size(mutpos{si}, 1);
      tot_num_mutations = tot_num_mutations + num_mutations(si);
      tmp = regexprep(M{si}.mut.newbase(midx),'^(.).*$','$1');
      newbase{si} = listmap(tmp, {'A', 'C', 'G', 'T'}); % 1=A 2=C 3=G 4=T 5=indel
      mut_categories{si} = M{si}.mut.categ_ignoring_null_categ(midx);
      mut_senses{si} = nan(num_mutations(si),1); % 1=synon  2=missense  3=nonsense/readthrough/splice  4=indel (can occur anywhere)
      idx = grepi('syn|silent',M{si}.mut.type(midx),1); mut_senses{si}(idx) = 1;
      idx = grepi('missense',M{si}.mut.type(midx),1); mut_senses{si}(idx) = 2;
      idx = grepi('nonsense|non.?stop|splice|read.?thr',M{si}.mut.type(midx),1); mut_senses{si}(idx) = 3;
      idx = grepi('ins|del',M{si}.mut.type(midx),1); mut_senses{si}(idx) = 4; newbase{si}(idx) = 5;
      if any(isnan(newbase{si})) || any(isnan(mut_senses{si}))
        fprintf('WARNING: nan in newbase or mut_senses: might cause error in gene %s\n',gname);
      end
    
      if P.mutsig2_use_sample_specific_coverage
        fprintf('Loading sample-specific coverage....\n');
        ssCov = nan(length(curr_genomic_array),length(midx));
        % find out what patients these mutations are in
        pidx = M{si}.mut.pat_idx(midx);
        % see if coverage is available for all those patients
        ok = nansub(M{si}.pat.fwb_ok,pidx);
        ssCov(:,ok) = get_from_fwb(M{si}.pat.fwb(pidx(ok)),chr,curr_genomic_array,M{si}.build,P.mutsig2_coverage_fwi);
        % for samples without coverage available, copy the median coverage from at least 100 samples
        if any(~ok)
          min_median_ct = P.min_num_samples_to_median_for_unknown_coverage;
          min_median_ct = min(min_median_ct,sum(M{si}.pat.fwb_ok));
          if sum(ok)>=min_median_ct      % if enough mutations, use the median coverage of the mutated samples
            medCov = ceil(median(ssCov(:,ok),2));
          else                           % otherwise supplement with randomly chosen patients
            nneed = min_median_ct-sum(ok);
            idx = setdiff(find(M{si}.pat.fwb_ok),pidx(ok));
            if nneed>length(idx), nneed=length(idx); end
            rp = randperm(length(idx));
            ii = idx(rp(1:nneed));
            otherCov = get_from_fwb(M{si}.pat.fwb(ii),chr,curr_genomic_array,M{si}.build,P.mutsig2_coverage_fwi);
            medCov = ceil(median([ssCov(:,ok) otherCov],2));
          end
          ssCov(:,~ok) = repmat(medCov,1,sum(~ok));
        end
        % impose coverage=true in radius around each mutation
        mask_radius = P.radius_to_impute_coverage_around_mutations;
        for i=1:length(midx)
          dnapos = curr_genomic_array(mutpos{si}(i));
          idx = find(abs(curr_genomic_array-dnapos)<=mask_radius);
          ssCov(idx,i) = 1;
        end
        % also compute summed coverage for display in generate.m
        coverage_track{si} = sum(ssCov,2);
        maxCov(si) = max(coverage_track{si});
        
      elseif P.mutsig2_use_power_calculation
        fprintf('Loading sample-specific coverage with raw read counts for power calculation....\n');
        P = impose_default_value(P, 'mutsig2_inwidth', 8);
        ssCov = nan(length(curr_genomic_array),length(midx));
        % find out what patients these mutations are in
        pidx = M{si}.mut.pat_idx(midx);
        pidx = unique(pidx);
        ssCov = nan(length(curr_genomic_array),length(pidx));
        % see if coverage is available for all those patients
        ok = nansub(M{si}.pat.power_ok,pidx);
        %%% For now use matrix of random numbers between 1 and 127, while files are copying 
        % if there is at least one sample with coverage available, extract coverage from those
        if any(ok)
          ssCov(:,ok) = get_from_fwb(M{si}.pat.power_files(pidx(ok)),chr,curr_genomic_array,M{si}.build,P.mutsig2_power_coverage_fwi);
        end 
        % for samples without coverage available, copy the median coverage from at least 100 samples
        if any(~ok)
          min_median_ct = P.min_num_samples_to_median_for_unknown_coverage;
          min_median_ct = min(min_median_ct,sum(M{si}.pat.power_ok));
          if sum(ok)>=min_median_ct      % if enough mutations, use the median coverage of the mutated samples
            medCov = ceil(median(ssCov(:,ok),2));
          else                           % otherwise supplement with randomly chosen patients
            nneed = min_median_ct-sum(ok);
            idx = setdiff(find(M{si}.pat.power_ok),pidx(ok));
            if nneed>length(idx), nneed=length(idx); end
            rp = randperm(length(idx));
            ii = idx(rp(1:nneed));
            otherCov = get_from_fwb(M{si}.pat.power_files(ii),chr,curr_genomic_array,M{si}.build,P.mutsig2_power_coverage_fwi);
            medCov = ceil(median([ssCov(:,ok) otherCov],2));
          end
          ssCov(:,~ok) = repmat(medCov,1,sum(~ok));
        end
        % for samples without coverage available, copy the median coverage from at least 100 samples        
        coverage_track{si} = sum(ssCov, 2) ~= 0;
        maxCov(si) = max(coverage_track{si});
      else   % use summed coverage profile
        
        if ~P.impute_full_coverage
          fileCov = org.broadinstitute.cga.tools.seq.FixedWidthBinary(M{si}.file.summed_cov_track);
          fileCov.setNullVal(0);
          coverage_track{si} = populate_coverage(curr_genomic_array, fileCov, chr);
          fileCov.close();
        else
          coverage_track{si} = ones(size(curr_genomic_array,2),1);
          if isfield(M{si},'np')
            coverage_track{si} = coverage_track{si} * M{si}.np; % (will be factored out later anyways)
          end
        end
        if sum(coverage_track{si}) == 0
          fprintf('No coverage found in gene. Skipping it...\n');
          fprintf('[REPORT]  %-11s   SKIP: no coverage\n',gname);
          skip_this_gene = true;
          break;
        end

        % impose summed_coverage=at_least_median in radius around each mutation
        medcov = median(coverage_track{si}(coverage_track{si}>0));
        mask_radius = P.radius_to_impute_coverage_around_mutations;
        for i=1:length(midx)
          dnapos = curr_genomic_array(mutpos{si}(i));
          idx = find(abs(curr_genomic_array-dnapos)<=mask_radius);
          coverage_track{si}(idx) = max(coverage_track{si}(idx),medcov);
        end
        
        %%  [ML 2011-09-20]
        %%  simplify coverage track
        %%  =======================
        %%  This step was necesasry to keep memory requirements under control.
        %%  (e.g. MUC16 in the pan-cancer 3000-exome set was generating a 15GB "throwable" matrix
        %%   and crashing the scatter-gather job)
        %%  To combat this problem, we will quantize coverage into a maximum number of bins.
        %%  (e.g. imposing limit of 100 bins reduced the MUC16 "throwable" matrix to 0.5GB)
        %%  Then finally, to keep the correct y-axis in the plots produced by generate(),
        %%  we tell it a scaling factor "coverage_track_factor" to use.
        
        maxCov(si) = max(coverage_track{si});
        if maxCov(si) > P.mutsig2_max_coverage_bins
          coverage_track_factor(si) = maxCov(si) / P.mutsig2_max_coverage_bins;
          coverage_track{si} = round(coverage_track{si} / coverage_track_factor(si));
          maxCov(si) = max(coverage_track{si});
        end
      end

    end % next si

    if tot_num_mutations < 2
        fprintf('Fewer than 2 mutations to analyze: skipping this gene...\n');
        fprintf('[REPORT]  %-11s   SKIP: fewer than two mutations\n',gname);
        skip_this_gene = true;
    end
    
    if skip_this_gene
      continue
    end

    % use whichever (hg18/hg19) track exists (first one, if both exist)
    conservation_to_use = conservation_track{bi(1)};
    context_and_sense_to_use = context_and_sense_track{bi(1)};
    context_and_sense_categs_to_use = context_and_sense_categs{bi(1)};

    % GENERATE LIST OF THROWABLE MUTATIONS

    fprintf('Generating list of throwable mutations...\n');
    throwable = cell(nsets,1);
    error_found = false;
    for si=1:nsets

        %%% modified 2012-01-17
        %%% -- when permuting, preserve:   category (from mutation context, as we already have been doing)
        %%%                                sense (syn/mis/nons+splice)--> but indels can go anywhere

        Q = collapse_context_and_sense_categories(context_and_sense_categs_to_use,M{si}.categ);

        if ~P.mutsig2_use_sample_specific_coverage  && ~P.mutsig2_use_power_calculation % permutations can cross patient boundaries
          [ucs tmp throw_flavor{si}] = unique([mut_categories{si} mut_senses{si}],'rows');
          nflavors = size(ucs,1);
          throwable{si} = cell(nflavors,1);
          for fi=1:nflavors
            c = ucs(fi,1);
            s = ucs(fi,2);
            if s<4
              [context_and_sense_idx newbase_idx] = find(Q(:,c,:)==s);
            else % indel: doesn't need to match syn/mis/nons
              [context_and_sense_idx newbase_idx] = find(Q(:,c,:)>0);
            end
            throwable{si}{fi} = populate_cont_pairs(coverage_track{si}, context_and_sense_to_use, ...
                                                     [context_and_sense_idx newbase_idx]);
          end
        elseif P.mutsig2_use_power_calculation
          nflavors = length(midx);
          throw_flavor{si} = (1:length(midx))';
          throwable{si} = cell(nflavors,1);
          for mi=1:length(midx)
            c = mut_categories{si}(mi);
            s = mut_senses{si}(mi);
            pix = find(pidx == M{si}.mut.pat_idx(midx(mi)));
            if s<4
              [context_and_sense_idx newbase_idx] = find(Q(:,c,:)==s);
            else   % indel: doesn't need to match syn/mis/nons
              [context_and_sense_idx newbase_idx] = find(Q(:,c,:)>0);
            end
            csii = listmap(context_and_sense_to_use,context_and_sense_idx);
            elig = find(~isnan(csii) & coverage_track{si});
            throwable{si}{mi} = [elig newbase_idx(csii(elig))];            
            try 
              tot_count = ssCov(find(M{si}.mut.start(midx(mi)) == curr_genomic_array), pix);
              alt_count = floor(tot_count*M{si}.mut.i_tumor_f(midx(mi)));
              ref_count = tot_count-alt_count;
              throwable{si}{mi} = power_calculation(throwable{si}{mi}, ssCov(:, pix), alt_count, ref_count, P);
            catch 
              keyboard 
            end 
            
          end
        else  % using sample-specific coverage, i.e. permutation within patients
          nflavors = length(midx);
          throw_flavor{si} = (1:length(midx))';
          throwable{si} = cell(nflavors,1);
          for mi=1:length(midx)
            c = mut_categories{si}(mi);
            s = mut_senses{si}(mi);
            if s<4
              [context_and_sense_idx newbase_idx] = find(Q(:,c,:)==s);
            else   % indel: doesn't need to match syn/mis/nons
              [context_and_sense_idx newbase_idx] = find(Q(:,c,:)>0);
            end
            csii = listmap(context_and_sense_to_use,context_and_sense_idx);
            elig = find(~isnan(csii) & coverage_track{si});
            throwable{si}{mi} = [elig newbase_idx(csii(elig))];
          end
        end
        for fi=1:nflavors
          if isempty(throwable{si}{fi})
            fprintf('ERROR: failed to generate throwable pairs for flavor %d in gene %s\n',fi,gname);
            error_found = true;
          end
        end
        if ~isempty(mutpos{si})
          if any(mutpos{si}<1 | isnan(mutpos{si}) | mutpos{si}>length(context_and_sense_to_use) | ...
                 mutpos{si}>length(conservation_to_use) | mutpos{si}>length(coverage_track{si}))
            fprintf('ERROR!\n');
            error_found = true;
            fprintf('[REPORT]  %-11s   SKIP: error in gene\n',gname);
            break  % skip this gene
          end
        end
    end
    
    if error_found
      fprintf('ERROR FOUND\n');
      continue  % skip this gene
    end

    % PERFORM PERMUTATIONS
    fprintf('Starting permutations...\n');

    try

      [p_clust, p_cons, p_joint, nperm, eff_clust, eff_cons] = generate(...
          genelength, nsets, P.setnames, tot_num_mutations, conservation_to_use, throw_flavor, mutpos, ...
          chr, gname, newbase, coverage_track, coverage_track_factor, throwable, maxCov, P);

    catch me
      fprintf('ERROR in permutations for %s\n',gname);
      for i=1:length(me.stack),disp(me.stack(i));end
      disp(me.message);
      disp(me);
      continue
    end

    R.nperm(ridx) = nperm;
    R.p_clust(ridx) = p_clust;
    R.p_cons(ridx) = p_cons;
    R.p_joint(ridx) = p_joint;
    R.eff_clust(ridx) = eff_clust;
    R.eff_cons(ridx) = eff_cons;
    look(R,ridx)
    ridx=ridx+1;

end

end
