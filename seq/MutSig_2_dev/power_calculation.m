function new_throwable = power_calculation(curr_throwable, patCov, alt_count, ref_count, P) 
% For a current mutation, the function takes the throwable sites, the 
% coverage for its respective patient, the allele counts of the mutation, the method to use, 
% and a paramaters struct (P). It outputs for every throwable site the power to have detected the mutation there
% given its allelic fraction and the coverage at that site. 

P = impose_default_value(P, 'mutsig2_beta_binomial_k', 3);
P = impose_default_value(P, 'mutsig2_lod_c', 0);
P = impose_default_value(P, 'mutsig2_lod_e', 0.001);
P = impose_default_value(P, 'mutsig2_lod_thresh', 6.3);
P = impose_default_value(P, 'smooth', 0);
P = impose_default_value(P, 'mutsig2_power_calculation_method', 'lod_score');
method = P.mutsig2_power_calculation_method;

%% Calculate power depending on method 
cov_counts = patCov(curr_throwable(:, 1));
if strcmp(method, 'beta_binomial')
  k = P.mutsig2_beta_binomial_k;
  power = beta_binomial_power(ref_count, alt_count, cov_counts, k);
elseif strcmp(method, 'lod_score') 
  c = P.mutsig2_lod_c; 
  e = P.mutsig2_lod_e; 
  depth = cov_counts;
  delta = alt_count/(ref_count + alt_count);
  lod_thresh = P.mutsig2_lod_thresh; 
  smooth = P.smooth;
  power = arrayfun(@calc_tumor_power, depth, repmat(e, length(depth), 1), repmat(lod_thresh, length(depth), 1), repmat(delta, length(depth), 1), repmat(c, length(depth), 1), repmat(smooth, length(depth), 1));
  if ~P.mutsig2_remove_nan_allelic_frac
    power(isnan(power)) = 1;
  end

end 

% Repopulate throwable matrix based on calculated power at each throwable site
try 
  new_throwable = nan(sum(floor(power*10)), 2); 
catch 
  keyboard 
end 
ix = 1; 
power(power < 0) = 0;
for i = 1:size(curr_throwable, 1)
  disp(i)
  try 
    new_throwable(ix:ix+floor(power(i)*10)-1, :) = repmat(curr_throwable(i, :), floor(power(i)*10), 1);
    ix = ix + floor(power(i)*10);
  catch me
    keyboard 
  end
end 

%% calculate power matrix (for testing purposes only)
function power_mat = calculate_power_matrix(depths, deltas, qscore, lod, contam, smoothing) 
    power_mat = repmat(-1, length(depths), length(deltas));
    for i = 1:length(depths)
        for j = 1:length(deltas)
            power_mat(i,j) = calc_tumor_power(depths(i), 10^(-1*qscore/10), lod, deltas(j), contam, smoothing)
        end
    end


%% calculate log likelihood for lod power calculation
function likelihood = ll(n, alts, e, f)

a = (n-alts) * log10(f*e/3 + (1-f)*(1-e));
b = (alts) * log10(f*(1-e) + (1-f)*e/3);
likelihood = a+b;

%% calculate tumor lod 
function t_lod = tumor_lod(n, alts, e, c) 

f = alts / n;
t_lod = ll(n, alts, e, f) - ll(n, alts, e, min(f,c));

%% calculate the power using the lod calculation method 
function power = calc_tumor_power(depth, e, lod_thresh, delta, c, smoothing) 

power = repmat(0, length(depth), 1);
for j = 1:length(power) 

  n = depth(j);
  a = 0:n;
  if delta==0 
    p_alt_given_e_delta = delta*(1-e) + (1-delta)*e;
  else 
    p_alt_given_e_delta = delta*(1-e) + (1-delta)*e/3;
  end
  p = binopdf(a, n, p_alt_given_e_delta); 
%  lod = 0:n; 
  %  for i = 1:n+1
  %    lod(i) = tumor_lod(n, a(i), e, c);
  %  end  
  lod = arrayfun(@tumor_lod, repmat(n, n+1, 1), a', repmat(e, n+1, 1), repmat(c, n+1, 1))';
  pass = find(lod >= lod_thresh); 
  if length(pass) > 0 
    k = min(pass);
    x = 0;
    if smoothing
      x = 1 - (lod_thresh - lod(k - 1))/(lod(k) - lod(k-1)); 
    end
    power(j) = sum(p(pass)) + x*p(k-1); 
  else 
    power(j) = 0; 
  end 

end 


%% Calculate the power using beta binomial
function power = beta_binomial_power(ref_count, alt_count, cov_counts, n) 

p = nan(length(cov_counts), 1);
p = arrayfun(@beta_binomial_prob, repmat(0, length(cov_counts), 1), cov_counts, repmat(alt_count, length(cov_counts), 1), repmat(ref_count, length(cov_counts), 1));
for z = 1:n-1 
  p = p + arrayfun(@beta_binomial_prob, repmat(z, length(cov_counts), 1), cov_counts, repmat(alt_count, length(cov_counts), 1), repmat(ref_count, length(cov_counts), 1));
end 
power = 1 - p;

%% beta binomial probability 
function p = beta_binomial_prob(k,n,a,b) 

z = 1; 
if k>n
  p = 0; 
else 
 p = nchoosek(n, k)*beta(k+a+z,n-k+b+z)/beta(a+z,b+z);
end 
  