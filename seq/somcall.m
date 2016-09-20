function M = somcall(t_bamfile,n_bamfile,chr,start,stop,params)
% M = somcall(t_bamfile,n_bamfile,chr,start,stop)
%
% columns of M are:
%   (1) chr
%   (2) position (1-based)
%   (3) reference base
%   (4) variant base
%   (5) tumor LOD:  (ref/alt + alt/alt) vs. ref/ref
%   (6) normal LOD: ref/ref vs (ref/alt + alt/alt)
%
%  [filtering statistics:]
%   (7) fraction of overlapping reads that are mapping-quality-zero (in tumor)
%   (8)   "" (in normal)
%   (9) average number of mismatches in overlapping reads (in tumor)
%  (10)   "" (in normal)
%  (11) number of distinct weird-pair destinations among overlapping reads (in tumor)
%  (12)   "" (in normal)
%  (13) estimated probability that SNP is nonartifactual (0-1)
%
% print_mutation_list(M) to convert to human-readable table
%
% Gaddy Getz and Mike Lawrence 2009-08

if ~exist('params','var')
  params=[];
end
params=impose_default_value(params,'oldstyle',false);
params=impose_default_value(params,'quiet',0);
params=impose_default_value(params,'maxsz',2e6);
params=impose_default_value(params,'nparts',1);
params=impose_default_value(params,'method','classic');
params=impose_default_value(params,'classic_thresh',6.3);
params=impose_default_value(params,'pv_thresh',5e-7);
params=impose_default_value(params,'n_LOD_thresh',2.3);
params=impose_default_value(params,'K',[0.2 8 5]);  % filtering: inflection point
params=impose_default_value(params,'Z',[5 5 5]);    % filtering: sigmoidicity

% handle multiple parts
if params.nparts>1
  nparts=params.nparts;
  params=impose_default_value(params,'lsfdir',pwd);
  params=impose_default_value(params,'readlen',100);
  params=impose_default_value(params,'queue',[]);
  params=impose_default_value(params,'rusage',[]);
  
  sz=stop-start+2*params.readlen;
  csz=cumsum(sz);
  
  if (csz(end)-max(sz))/nparts > params.maxsz
    error('parts exceed max size');
  end
  
  p=cell(nparts,1);
  k=1;
  cursz=0;
  for i=1:nparts
    next_size=(csz(end)-cursz)/(nparts-i+1);
    idx=find(csz<=(next_size+cursz),1,'last');
    if isempty(idx) || idx<k
      error('can not allocate next part');
    end
    p{i}=k:idx;
    k=idx+1;
    cursz=csz(idx);
  end
  
  % p=get_parts(1:length(X.start),nparts);
  part_params=params;
  part_params.nparts=1;
  
  l=lsf(params.lsfdir);
  h=zeros(nparts,1);
  for i=1:nparts
    part_chr=chr(p{i});
    part_start=start(p{i});
    part_stop=stop(p{i});
    [l,h(i)]=bsub(l,{'M'},'somcall',{t_bamfile,n_bamfile,part_chr,part_start,part_stop,part_params},params.queue,params.rusage);
  end
  [l,res]=wait(l); % wait for all
  
  Ms=cell(nparts,1);
  for i=1:nparts
    Ms({i})=res{h(i)}.M;
  end
  M=cat(1,Ms{:});
  return
end

% handle multiple chromosomes
if length(chr)>1
  [uchr,uchri,uchrj]=unique(chr);
  if length(uchr)>1
    for i=1:length(uchr)
      idx=find(uchrj==i);
      Ms{i}=somcall(t_bamfile,n_bamfile,uchr(i),start(idx),stop(idx),params);
    end
    M=cat(1,Ms{:});
    return
  else
    chr=uchr;
  end
end

subfprintf('SOMCALL starting at %s\n',datestr(now));
tt0=tic;

%              ref/ref  ref/alt alt1/alt1 alt1/alt2 hom_N het_N
n_priors=log10([ 0.25     0.25    0.25      0.25      1/2   1/2   ]); % no priors 
t_priors=log10([ 0.25     0.25    0.25      0.25      1/2   1/2   ]); % no priors 

n_params = params;
t_params = params;

n_params.filter_mapq0 = 0;
t_params.filter_mapq0 = 1;

t_params.trim_offtarget_bases = false;    % time savings is less than the time it takes to trim
n_params.trim_offtarget_bases = true;    % time+memory savings is huge

switch params.method
 case 'pv';
  t_params.calc_mu_var=true;
  n_params.calc_mu_var=false;
 case 'classic'
  t_params.calc_mu_var=false;
  n_params.calc_mu_var=false;
 otherwise
  error('no such method');
end

subfprintf('\nCalc_models for T:');
[t_LL,t_refref_model,t_ref,t_coverage,t_R,t_B,t_ri,t_bi,t_refi,t_LLmu,t_LLvar]=calc_models(t_bamfile,chr,start, ...
                                                  stop,t_priors,t_params);
t_ref_N_idx=find(t_ref>4);

if isempty(t_R) | isempty(t_B) 
  subfprintf('No tumor data!\n');
  M = zeros(0,13);
  return;
end

tt1=toc(tt0);

% Tumor:  LOD (ref/alt + alt/alt)vs ref/ref  for all alts
subfprintf('Collapse models for T...');
[t_LL_model,t_refref_LL]=calc_collapsed_models(t_LL,t_refref_model,t_ref);

[mx_t_LL_model,mx_t_LL_model_i]=max(t_LL_model,[],1);
t_LOD=mx_t_LL_model-t_refref_LL;
t_LOD(t_ref_N_idx)=-Inf;

switch params.method
 case 'pv'
  t_LL_cutoff=t_LLmu+sqrt(t_LLvar)*norminv(params.pv_thresh,0,1);
  candidate_idx=find(t_refref_LL<t_LL_cutoff');
 case 'classic'
  candidate_idx=find(t_LOD>params.classic_thresh);
 otherwise
  error('no such method');
end
  

if isempty(candidate_idx)
  subfprintf('\n\nNo candidates!\n');
  M = zeros(0,13);
  return;
end

subfprintf('\n\n%d candidates\n',length(candidate_idx));

[cand_region,cand_pos]=map_to_intervals(candidate_idx,start,t_refi);

tt2=toc(tt0);

subfprintf('\nCalc_models for N:');
[n_LL,n_refref_model,n_ref,n_coverage,n_R,n_B,n_ri,n_bi,n_refi]=calc_models(n_bamfile,chr,cand_pos,cand_pos,n_priors,n_params);
n_ref_N_idx=find(n_ref>4);

tt3=toc(tt0);

% Normal: LOD ref/ref vs (ref/alt + alt/alt)  for alt in tumor
subfprintf('Collapse models for N...'); 
[n_LL_model,n_refref_LL]=calc_collapsed_models(n_LL,n_refref_model,n_ref);
n_LOD=n_refref_LL-n_LL_model(mx_t_LL_model_i(candidate_idx)+4*(0:(length(n_ref)-1)));
n_LOD(n_ref_N_idx)=-Inf;

ref_in_norm_idx=find(n_LOD>params.n_LOD_thresh);

[somatic_region,somatic_pos]=map_to_intervals(ref_in_norm_idx,cand_pos,n_refi);

somatic=candidate_idx(somatic_region);

subfprintf('\n\n%d somatic mutations\n',length(somatic));

tt4=toc(tt0);

% FILTERING: STATISTICS

subfprintf('\nFiltering...'); 

if ~isempty(somatic)
  havebooms = ~isempty(grep('\.boom$',t_bamfile,1)) & ~isempty(grep('\.boom$',n_bamfile,1));
  if havebooms   % faster method (pre-computed stats in booms)
    t_stem = [t_bamfile '/chr' num2str(chr) '.'];
    n_stem = [n_bamfile '/chr' num2str(chr) '.'];
    t_fmq0 = get_block([t_stem 'fmapqz'],'byte',somatic_pos-1)/100;
    n_fmq0 = get_block([n_stem 'fmapqz'],'byte',somatic_pos-1)/100;
    t_avgnmm = get_block([t_stem 'avgnmm'],'byte',somatic_pos-1)/10;
    n_avgnmm = get_block([n_stem 'avgnmm'],'byte',somatic_pos-1)/10;
    t_nuwp = get_block([t_stem 'nuwp'],'byte',somatic_pos-1);
    n_nuwp = get_block([n_stem 'nuwp'],'byte',somatic_pos-1);
  else    % old, slower method (for bams)
    z = nan(length(somatic),1);
    t_fmq0 = z; n_fmq0 = z;
    t_avgnmm = z; n_avgnmm = z;
    t_nuwp = z; n_nuwp = z;
    n_ri2 = [0;n_ri];
    t_ri2 = [0;t_ri];
    for i=1:length(somatic)
      tr1 = t_ri2(cand_region(somatic_region(i)))+1; tr2 = t_ri2(cand_region(somatic_region(i))+1); 
      nr1 = n_ri2(somatic_region(i))+1; nr2 = n_ri2(somatic_region(i)+1);
      tidx = (tr1-1)+find(t_R(tr1:tr2,4)<=somatic_pos(i) & t_R(tr1:tr2,5)>=somatic_pos(i));
      nidx = (nr1-1)+find(n_R(nr1:nr2,4)<=somatic_pos(i) & n_R(nr1:nr2,5)>=somatic_pos(i)); 
      t_fmq0(i) = mean(t_R(tidx,8)==0);  % fraction MAPQZ=0 at this position
      n_fmq0(i) = mean(n_R(nidx,8)==0);
      t_avgnmm(i) = mean(t_R(tidx,7));   % average number of mismatches
      n_avgnmm(i) = mean(n_R(nidx,7));
      t_nuwp(i) = count_uwp(t_R,tidx,chr,somatic_pos(i));  % number of unique weird-pairs destinations
      n_nuwp(i) = count_uwp(n_R,nidx,chr,somatic_pos(i));
    end
  end
else
    t_fmq0 = [];
    n_fmq0 = [];
    t_avgnmm = [];
    n_avgnmm = [];
    t_nuwp = [];
    n_nuwp = [];
end


% FILTERING: SCORE

% tests=
%    (1)   fraction of overlapping reads with MAPQ=0
%    (2)   average number of mismatches in overlapping reads
%    (3)   number of unique weird-pair destinations among overlapping reads

% scoring=
%    (1) calculate probability that the SNP is real, based on each test
%        -- from sigmoidal curve with two test-specific parameters
%    (2) multiply probabilities from each test

if ~isempty(somatic)
  t_S = [t_fmq0 t_avgnmm t_nuwp];
  n_S = [n_fmq0 n_avgnmm n_nuwp];
  tmp =  bsxfun(@power,bsxfun(@rdivide,t_S,params.K),params.Z);
  t_P = 1 - (tmp ./ (1+tmp));
  tmp =  bsxfun(@power,bsxfun(@rdivide,n_S,params.K),params.Z);
  n_P = 1 - (tmp ./ (1+tmp));
  P = [t_P n_P];
  filter = prod(P,2);
else
  filter = [];
end

M = [...
    repmat(chr,length(somatic),1),...
    somatic_pos,...
    t_ref(somatic)',...
    mx_t_LL_model_i(somatic)',...
    t_LOD(somatic)',n_LOD(somatic_region)',...
    t_fmq0,n_fmq0,...
    t_avgnmm,n_avgnmm,...
    t_nuwp,n_nuwp,...
    filter...
];

tt5=toc(tt0);

subfprintf('  %d/%d have filter>0.5\n',sum(filter>0.5),length(somatic));

subfprintf('\nSOMCALL finished at %s\n',datestr(now));
subfprintf(['calcmodT %.2f   collapseT %.2f   calcmodN %.2f   collapseN %.2f   '...
    'filter %.2f   total %.2f s\n'],tt1,tt2-tt1,tt3-tt2,tt4-tt3,tt5-tt4,tt5);

return



% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
%  function nuwp = subfunction_nuwp(R,idx,pos)
%    wp = R(idx,10:12);                                    % chr start strand
%    wp = wp(wp(:,1)>0,:);                                 % remove unpaired and unmapped-pair
%    nw = find(wp(:,1)==chr & abs(wp(:,2)-pos)<10000);     % remove non-weird
%    wp = wp(setdiff(1:size(wp,1),nw),:);
%    wp(:,2) = round(wp(:,2) / 2000);                      % smooth to windows of 2kB
%    nuwp = size(unique(wp,'rows'),1);
%  end

  function subfprintf(str,varargin)
    if ~(params.quiet), fprintf(str,varargin{:}); end
  end     

end % of main function
