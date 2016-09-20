function nnon_pass = sns(nnon,nsil,fnon_exp)

if ~exist('fnon_exp','var'), fnon_exp=3/4; end

[fnon_obs fnon_obs_ci] = binofit(nnon,nnon+nsil);
fnon_obs_stdev = max((fnon_obs-fnon_obs_ci(1)),(fnon_obs_ci(2)-fnon_obs))/1.96;

% given observed number of silent mutations and fnon_exp,
%    expected number of nonsilent

ns_s_ratio_exp = 1 / ((1/fnon_exp)-1);
nnon_exp = nsil * ns_s_ratio_exp;



alpha = 0.05;
for n=0:nsil+10
  p = binopdf(n,n+nsil,fnon_exp);
%  [m ci] = binofit(n,n+nsil,fnon_exp);
  fprintf('%.0f non + %.0f sil    %0.2f\n',n,nsil,p);
end

return


alpha = 0.5;
for nnon_pass=0:nnon
  [m ci] = binofit(nnon-nnon_pass,nsil+nnon-nnon_pass,fnon_exp);
  fprintf('remove %.0f    %.0f non + %.0f sil    %0.2f - %0.2f\n',nnon_pass,nnon-nnon_pass,nsil,ci(1),ci(2));
end



