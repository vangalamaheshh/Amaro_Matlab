function [p_clust, p_cons, p_joint, nperm, eff_clust, eff_cons] = generate(...
    genelength, nsets, setnames, tot_nmuts, conservation_track, throw_flavor, mutpos, ...
    chr, gname, newbase, coverage_track, coverage_track_factor, throwable, maxCov, P)
%For each of the maxperm permutations:
%   1. Throws nmuts random mutations and places them in an array
%   2. Calculates the relative distances of all of the random mutations and
%       and uses them to calculate a metric.
%   3. Proceeds to populate an array of size maxperm with the obtained
%       metrics.

if ~exist('P','var'), P=[]; end
P = impose_default_value(P, 'mutsig2_maxperm', 5000000);  % maximum permutations per gene
P = impose_default_value(P, 'mutsig2_maxsec', inf);       % maximum seconds per gene
P = impose_default_value(P, 'mutsig2_theta', 0.1);        % keep going until reach this fractional confidence in p-value
P = impose_default_value(P, 'mutsig2_imagedir', []);      % if unspecified, then won't save images to disk
P = impose_default_value(P, 'mutsig2_repress_image_output', true);
P = impose_default_value(P, 'mutsig2_simulation1',false);
P = impose_default_value(P, 'mutsig2_keyboard_every_report',false);
P = impose_default_value(P, 'mutsig2_clustering_metric',1);  % 1=K-S
P = impose_default_value(P, 'mutsig2_convert_conservation_values_to_ranks',false);

if isfield(P,'mutsig2_uncertainty_tolerance')
  if isfield(P,'mutsig2_min_effect_size'), error('please use only one'); end
  fprintf('*** NOTE  --> P.mutsig2_uncertainty_tolerance is deprecated.  P.use mutsig2_min_effect_size instead.\n');
  fprintf('          --> replacing P.mutsig2_uncertainty_tolerance=%d to P.use mutsig2_min_effect_size instead=%d\n',...
          P.mutsig2_uncertainty_tolerance, 1+P.mutsig2_uncertainty_tolerance);
end
P = impose_default_value(P, 'mutsig2_min_effect_size',1.01);
if P.mutsig2_min_effect_size<1, error('P.mutsig2_min_effect_size must 1.00 or greater'); end

if ~isempty(P.mutsig2_imagedir), ensure_dir_exists(P.mutsig2_imagedir); end

DRAW_PIC_ONLY_AT_END = false;
if P.mutsig2_repress_image_output, DRAW_PIC_ONLY_AT_END = true; end

if tot_nmuts < 2
  fprintf('Need at least two mutations.\n');
  return
end
obs_mutpos = cat(1,mutpos{:});

% BUFFERS FOR PERMUTATION CREATION
perm_mutpos = cell(nsets,1);
nflavor = nan(nsets,1);
flavor_counts = cell(nsets,1);
one_per_flavor = false(nsets,1);
nthrowable = cell(nsets,1);
for si=1:nsets
  nflavor(si) = length(throwable{si});
  flavor_counts{si} = histc(throw_flavor{si},1:nflavor(si));
  nthrowable{si} = cellfun('length',throwable{si});
  one_per_flavor(si) = all(flavor_counts{si}==1);
end

%%%%%%%%%% RUN SIMULATION FIRST? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if P.mutsig2_simulation1
  sub_run_simulation();  % (at end of this file)
  % then proceed as normal, to return the actual p-values
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RANGE OF METRICS

if P.mutsig2_convert_conservation_values_to_ranks
  [tmp ord] = sort(conservation_track);
  [tmp conservation_track] = sort(ord);
end

min_clust = 0;
max_clust = 1;
min_cons = min(conservation_track); 
max_cons = max(conservation_track);

% bins for joint distribution
nbins = 100;
binsize_clust = (max_clust-min_clust)/(nbins-1);
binsize_cons = (max_cons-min_cons)/(nbins-1);

% CALCULATE METRICS FOR THE OBSERVED MUTATIONS

% observed conservation
obs_cons_unreduced = nanmean(conservation_track(obs_mutpos));
obs_cons = obs_cons_unreduced / P.mutsig2_min_effect_size;
obs_cons_bin = 1+floor((obs_cons - min_cons)/binsize_cons);

% observed clustering
if P.mutsig2_clustering_metric==1
  % Kolmogorov-Smirnoff test
  obs_curve = fast_calculate_curve(obs_mutpos,genelength);
  mean_curve = curve_mean(genelength);
  obs_clust_unreduced = max(obs_curve - mean_curve);
else
  obs_clust_unreduced = new_clustering_statistic(obs_mutpos,genelength,P.mutsig2_clustering_metric,true);  
end
obs_clust = obs_clust_unreduced / P.mutsig2_min_effect_size;
obs_clust_bin = 1+floor((obs_clust - min_clust)/binsize_clust);

% observed joint score
obs_joint = obs_clust + obs_cons; % (temporarily)

% BUFFERS FOR PERMUTATION RESULTS
perm_clust = zeros(P.mutsig2_maxperm,1);
perm_cons = zeros(P.mutsig2_maxperm,1);
k_clust = 0;
k_cons = 0;
k_joint = 0;
joint_hist = zeros(nbins,nbins);

% PERMUTATIONS

tictemp = tic;

nperm = 0;

first_check = 1000;
check_every = 10000;

tic1=0;
tic2=0;
tic3=0;
tic4=0;

aa=tic;

finished = false;
while(~finished)

  %%%%%%%%%%%%%%%%%%%%%%%%  DO A PERMUTATION

  nperm = nperm + 1;

  tic4=tic4+toc(aa);

  aa=tic;
  % RANDOMLY THROW MUTATIONS
  sub_generate_perm_mutpos();
  all_perm_mutpos = cat(1,perm_mutpos{:}); % collapse across sets
  tic1=tic1+toc(aa);

  aa=tic;
  % CALCULATE CONSERVATION METRIC FOR THIS PERMUTATION
  perm_cons(nperm) = nanmean(conservation_track(all_perm_mutpos(:,1)));
  tic2=tic2+toc(aa);

  aa=tic;
  % CALCULATE CLUSTERING METRIC FOR THIS PERMUTATION
  if P.mutsig2_clustering_metric==1
    % Kolmogorov-Smirnoff test
    perm_curve = fast_calculate_curve(all_perm_mutpos(:,1), genelength);
    perm_clust(nperm) = max(perm_curve - mean_curve);
  else
    perm_clust(nperm) = new_clustering_statistic(all_perm_mutpos(:,1),genelength,P.mutsig2_clustering_metric);
  end

  % INCREMENT HISTOGRAMS
  if perm_cons(nperm)>=obs_cons, k_cons=k_cons+1; end
  if perm_clust(nperm)>=obs_clust, k_clust=k_clust+1; end
  bin_cons = 1+floor((perm_cons(nperm) - min_cons)/binsize_cons);
  bin_clust = 1+floor((perm_clust(nperm) - min_clust)/binsize_clust);
  joint_hist(bin_clust,bin_cons) = joint_hist(bin_clust,bin_cons)+1;

  tic3=tic3+toc(aa);

  aa=tic;

  %%%%%%%%%%%%%%%%%%%%%%%%  CHECK P-VALUES?

  if nperm~=first_check && mod(nperm,check_every)>0 && nperm<P.mutsig2_maxperm
    continue   % not time to check p-values
  end

  % CALCULATE MARGINAL P-VALUES and effect sizes
  [p_clust ci_ratio_clust] = calc_pval_and_ci_ratio(k_clust,nperm);
  eff_clust = obs_clust_unreduced / mean(perm_clust(1:nperm));
  [p_cons ci_ratio_cons] = calc_pval_and_ci_ratio(k_cons,nperm);
  eff_cons = (max_cons-mean(perm_cons(1:nperm))) / (max_cons-obs_cons_unreduced);

  % CALCULATE JOINT P-VALUE
  landfill_dist = landfill(joint_hist/nperm);
  bin_score = -log10(landfill_dist);
  obs_score = bin_score(obs_clust_bin,obs_cons_bin);
  k_joint = sum(joint_hist(bin_score>=obs_score));
  [p_joint ci_ratio_joint] = calc_pval_and_ci_ratio(k_joint,nperm);
    
  %%%%%%%%%%%%%%%%%%%%%%%%  UPDATE TEXT REPORT

  tt = toc(tictemp);
  report_type = 2;
  if report_type==1
    fprintf('[REPORT]  %-11s nmuts %-4d len %-5d kcons %-5d kclust %-5d kjoint %-5d nperm %-7d %4.1f sec  (%.0f perm/sec)',...
            gname,tot_nmuts,genelength,k_cons,k_clust,k_joint,nperm,tt,nperm/tt);
    ttot=(tic1+tic2+tic3+tic4);
    ft1=tic1/ttot;ft2=tic2/ttot;ft3=tic3/ttot;ft4=tic4/ttot;
    fprintf(' %.0f/%.0f/%.0f/%.0f',ft1*100,ft2*100,ft3*100,ft4*100);
    fprintf('\n');
  elseif report_type==2
    fprintf(['[REPORT]  %-11s nmuts %-4d len %-5d nperm %-7d  (%0.f perm/sec)'...
              '  CONS k %-5d p %-0.6f eff %-2.2f'...
             '  CLUST k %-5d p %-0.6f eff %-2.2f'...
             '  JOINT k %-5d p %-0.6f\n'],...
            gname,tot_nmuts,genelength,nperm,nperm/tt,...
            k_cons,p_cons,eff_cons,...
            k_clust,p_clust,eff_clust,...
            k_joint,p_joint);
  end
keyboard
  if P.mutsig2_keyboard_every_report, keyboard; end

  %%%%%%%%%%%%%%%%%%%%%%%%  CHECK STOPPING CONDITIONS

  % CHECK IF WE HAVE CONVERGED ON AN ANSWER
  max_ci_ratio = max([ci_ratio_joint,ci_ratio_cons,ci_ratio_clust]);
  if (max_ci_ratio<=P.mutsig2_theta), finished = true; end

  % CHECK TIME ELAPSED
  if (tt>=P.mutsig2_maxsec), fprintf('Max time elapsed.\n'); finished = true; end
  
  % CHECK MAXPERM
  if (nperm>=P.mutsig2_maxperm), fprintf('Max permutations reached.\n'); finished = true; end

  if ~finished && DRAW_PIC_ONLY_AT_END
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%  DRAW/UPDATE PLOT
  clf
  set(gcf,'position',[10 10 1420 780]);
  gene_name = strcat('Name: ', gname, '  nperm = ',num2str(nperm));
  chr_number = strcat('Chr Number: ', num2str(chr));
  mutation_name = strcat('Mutations: ', num2str(obs_mutpos'));
  p_mle_text = strcat('Clustering P value: ', num2str(p_clust));
  p_cons_text = strcat('Conservation P value: ', num2str(p_cons));
  p_joint_text = strcat('Joint P value: ', num2str(p_joint),'  max_ci_ratio = ',num2str(max_ci_ratio));
  
  subplot(2,4,1)
  hold on
  cccs = distinct_colors(nsets);
  for si=1:nsets
    ccc = cccs(si,:);
    plot(coverage_track{si},'color',ccc);
    mty = -inf;
    for m = 1:size(mutpos{si},1)
      ty = maxCov(si)+m/10;
      text(mutpos{si}(m), ty, 'O','color',ccc)
      mty = max(mty,ty);
    end
  end
        
  xlabel('Length of Gene');
  ylabel('Coverage');
  title('Coverage Track');
  if nsets>1 && length(unique(setnames))>1
    if nsets>2, loc='eastoutside'; else loc='southwest'; end
    legend(setnames,'location',loc);
  end
  yl = ylim;
  yl(2) = max(yl(2),mty);
  yl(1) = min(yl(1),yl(2)/2);
  ylim(yl);
  xlim([1 length(conservation_track)]);
  hold off
    
  subplot(2,4,5);
  scatter(1:length(conservation_track),conservation_track,5);
  xlabel('Length of Gene');
  ylabel('Conservation Value');
  xlim([1 length(conservation_track)]);
  title('Conservation Track for this Gene');
    
  subplot(2,4,2);
  hist(perm_clust(1:nperm),min_clust:binsize_clust:max_clust);
  mx = max(histc(perm_clust(1:nperm), min_clust:binsize_clust:max_clust));
  line([obs_clust obs_clust], [0 mx], 'Color', [1.0 0.0 0.0]);
  text(obs_clust, mx*0.9,num2str(p_clust));
  xlim1 = min_clust;
  xlim2 = max(max_clust,obs_clust);
  xlim([xlim1 xlim2]);
  xlabel('Clustering Score');
  ylabel('Number of Permutations');
  title('Marginal Distribution of Clustering Scores')
        
  subplot(2,4,7);
  hist(perm_cons(1:nperm),min_cons:binsize_cons:max_cons);
  line([obs_cons obs_cons_bin], [0 mx], 'Color', [1.0 0.0 0.0]);
  text(obs_cons, mx*0.9, num2str(p_cons));
  xlim([min(conservation_track),max(conservation_track)]);
  xlabel('Conservation Score');
  ylabel('Number of Permutations');
  title('Marginal Distribution of Conservation Scores');
      
  subplot(2,4,6);
  hist(conservation_track,50);
  xlim([min(conservation_track),max(conservation_track)]);
  xlabel('Conservation Value');
  ylabel('Number of Conservation Scores in Track');
  title('Conservation Distribution For This Gene');
  
  subplot(2,4,3);
  imagesc(joint_hist);
  colorbar;
  line(xlim, [obs_clust_bin obs_clust_bin], 'Color', [1.0 0.0 0.0])
  line([obs_cons_bin obs_cons_bin], ylim, 'Color', [1.0 0.0 0.0])
  xlabel('Conservation Score');
  ylabel('Clustering Score');
  title('Joint Distribution of Clustering and Conservation Scores');
  
  subplot(2,4,4);
  cap=8;
  imagesc(min(cap,landfill_dist));
  colorbar;
  line(xlim, [obs_clust_bin obs_clust_bin], 'Color', [1.0 0.0 0.0])
  line([obs_cons_bin obs_cons_bin], ylim, 'Color', [1.0 0.0 0.0])
  xlabel('Conservation Score');
  ylabel('Clustering Score');
  title('Landfilled Joint Distribution');

  subplot(2,4,8);
%  plot(score_hist); (1D joint score histogram no longer explicitly calculated)
%  line([score_mapped_index score_mapped_index], [0 max(score_hist)], 'Color', [1.0 0.0 0.0]);
%  ylim([0 max(score_hist)*1.3]);
%  title('Log scores probability distribution');
  text(0,1,{gene_name, chr_number, mutation_name,...
                      p_mle_text, p_cons_text, p_joint_text},'interpreter','none');
  set(gca,'visible','off');

end   % next permutation

%%%%% ENDPOINT REACHED

% print plot to file (if P.mutsig2_imagedir was specified)
if ~isempty(P.mutsig2_imagedir) && ~P.mutsig2_repress_image_output
  set(gcf,'paperpositionmode','auto','papersize',[11 8.5])
  temp_name = strcat(gname, '.', num2str(tot_nmuts), '_', num2str(genelength));
  full_name = [P.mutsig2_imagedir,'/',temp_name, '_', num2str(nperm), '_', ...
               num2str(p_clust), '_', num2str(p_cons), '_', num2str(p_joint), '.png'];
  print_to_file(full_name);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sub_generate_perm_mutpos()

  for si=1:nsets
    if one_per_flavor(si) % slightly faster than old way, in cases where we need to throw one mutation of each "flavor"
      idx = ceil(nthrowable{si}.*rand(nflavor(si),1));
      perm_mutpos{si} = nan(length(idx),2);
      for j=1:length(idx)
        perm_mutpos{si}(j,:) = throwable{si}{j}(idx(j),:);
      end
    else   % old way
      tmp = cell(nflavor(si), 1);
      for j = 1:nflavor(si)
        tmp{j} = throwable{si}{j}(ceil(size(throwable{si}{j},1)*rand(flavor_counts{si}(j),1)),:);
      end
      perm_mutpos{si} = cat(1, tmp{:});
    end
  end

  % CHECKED: allocating perm_mutpos once and updating in place: was ~17% SLOWER!
  % ALSO CHECKED: randi() is slower than rand()

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sub_run_simulation()

  P = impose_default_value(P,'mutsig2_num_randomized_versions',1000);

  PP=P;
  PP.mutsig2_simulation1 = false;
  PP.mutsig2_repress_image_output = true;

  fprintf('MutSig2 simulation mode:\n');

  n = P.mutsig2_num_randomized_versions;
  pclust = nan(n,1);
  pcons = nan(n,1);
  pjoint = nan(n,1);

  for i=1:n

    fprintf('[%d/%d]  ', i, n);
    
    % STEP1: generate a randomized version using the same procedure as the permutations
    %        (need to alter "mutpos"=cDNA position, "newbase")

    sub_generate_perm_mutpos();
    for si=1:nsets
      Xmutpos{si} = perm_mutpos{si}(:,1);
      Xnewbase{si} = perm_mutpos{si}(:,2);
    end

    % STEP2: test this randomized version via permutations
  
    [pclust(i), pcons(i), pjoint(i)] = generate(...
        genelength, nsets, setnames, tot_nmuts, conservation_track, throw_flavor, Xmutpos,...
        chr, gname, Xnewbase, coverage_track, coverage_track_factor, throwable, maxCov, PP);
  
  end % next randomized version

  qq(pclust,pcons,pjoint);
  legend({'pclust','pcons','pjoint'},'location','northwest','fontsize',12);

  fprintf('Simulation complete.  Type "dbcont" to calculate actual numbers\n');

end % of sub_run_simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % of main function










