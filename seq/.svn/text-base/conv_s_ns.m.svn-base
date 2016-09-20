function [pval ratio_exp] = conv_s_ns(nnon,nsil,ns_s_ratio)

ncat = length(nnon);
if length(nsil)~=ncat, error('length(fNS~=ncat)'); end
if length(ns_s_ratio)~=ncat, error('length(fNS~=ncat)'); end

fNS = 1./(1+(1./ns_s_ratio));

ntot = nnon+nsil;

NS = [];
for c=1:ncat
  if ntot(c)>0
    if fNS(c)>=0 && fNS(c)<=1
      NS{end+1,1} = binopdf(0:ntot(c),ntot(c),fNS(c));
    else
      fprintf('Unexpected fNS = %f in analyze_silent_nonsilent_ratios\n',fNS);
      keyboard
    end
  end
end

% convolute to get final distribution
dist = [];
dist.prob = batch_convolute(NS);
dist.nnon = (0:slength(dist)-1)';
dist.nsil = slength(dist)-1-dist.nnon;
dist.ratio = dist.nnon ./ dist.nsil;

% MLE
[tmp idx] = max(dist.prob);
nnon_exp = dist.nnon(idx);
nsil_exp = dist.nsil(idx);
ratio_exp = dist.ratio(idx);

% p-value
nnon_obs_tot = sum(nnon);
nsil_obs_tot = sum(nsil);
ratio_obs = nnon_obs_tot / nsil_obs_tot;
idx = find(dist.ratio>=ratio_obs);

if isempty(idx)
  fprintf('Unexpected case in analyze_silent_nonsilent_ratios\n');
  keyboard
  pval = 0;
elseif idx==1
  pval = 1;
else
  pval = 1-sum(dist.prob(1:idx-1));
end
