function [m n pmum z] = ctreads(R,chr)

nbins = max(R(:,13));

% for each binned start position
%   find reads that start at that position with readstrand=0 and matestrand=1
%     copmute n = number of such reads
%     compute m = mean insert size
idx = find(R(:,10)==chr & R(:,6)==0 & R(:,12)==1 & R(:,11)>R(:,4) & R(:,11)-R(:,4)<10000);
m = nan(nbins,1);
n = zeros(nbins,1);
[u ui uj] = unique(R(idx,13));
for i=1:length(u), if ~mod(i,1000), fprintf('%d/%d ',i,length(u)); end
  jdx = idx(uj==i);
  n(u(i)) = length(jdx);
  m(u(i)) = mean(R(jdx,14));
end, fprintf('\n');

%  find reads that start at that position with readstrand=0 and unmapped pairmate
%      compute pmum = number of such reads
idx = find(R(:,6)==0 & R(:,10)==0);
pmum = zeros(nbins,1);
[u ui uj] = unique(R(idx,13));
for i=1:length(u), if ~mod(i,1000), fprintf('%d/%d ',i,length(u)); end
  jdx = idx(uj==i);
  pmum(u(i)) = length(jdx);
end, fprintf('\n');

%      compute z = number of mapping-quality-zero reads
idx = find(R(:,10)==chr & R(:,6)==0 & R(:,12)==1 & R(:,11)>R(:,4) & R(:,11)-R(:,4)<10000 & R(:,8)==0);
z = zeros(nbins,1);
[u ui uj] = unique(R(idx,13));
for i=1:length(u), if ~mod(i,1000), fprintf('%d/%d ',i,length(u)); end
  jdx = idx(uj==i);
  z(u(i)) = length(jdx);
end, fprintf('\n');
