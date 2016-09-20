function m = dRanger_calculate_expected_stdev_distrib(distrib,P)
% dRanger_calculate_expected_stdev_distrib(distrib,P)
%
% INPUT
%
% distrib = insert-size distribution
%
%   if distrib is a string, it should point to a *_isz directory or a *.isz file
%   containing the empirical distribution
% 
%   if distrib is numeric, it should be a matrix characterizing
%   the distribution as a sum of Gaussians:
%
%          each row:  frac mean std
%          e.g.       0.5  250  15
%                     0.5  400  20
%
%   if distrib is a struct, it should have the following fields:
%          .dat = matrix of insert sizes (rows = isz, columns = lanes)
%          .mx  = 
%
% P.maxstdev = max stdev to calc P for
%          default = 200
%
% P.setsizes = set of possible sizes of the set of pairs
%          default = [2:6 8 10 12 15 20 30]
%
% P.nperms = number of random permutations to take:
%          1/nperms = lower bound on p-values that can be calculated
%          default = 1e7
%
% P.readlen = read length
%          default = 101
%
% P.maxpeek = maximum number of basepairs by which a read is allowed
%                   to peek across the breakpoint and still be aligned
%          default = 8
%
% OUTPUT
%
% m
%     number of rows = length(P.setsizes)
%     each row:
%         <setsize>   <mean_of_std>   <std_of_std>    <mean_of_range>    <std_of_range>
%
%     range=(max-min)
%     
% Mike Lawrence 2009-10

if ~exist('distrib','var'), error('Required: distrib'); end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'maxstdev', 200);
P = impose_default_value(P,'setsizes', [2:6 8 10 12 15 20 30]);
P = impose_default_value(P,'nperms', 1e6);
P = impose_default_value(P,'readlen', 101);
P = impose_default_value(P,'maxpeek',10);

fprintf('dRanger_calculate_expected_stdev_distrib\n');

nrows = 1e8;
setsizes = P.setsizes;

if isnumeric(distrib)

  % distribution expressed as sum of Gaussians
  if sum(distrib(:,1))~=1, error('frac must sum to 1'); end

  % choose cen = max(mean+6std)
  range = distrib(:,2) + 6*distrib(:,3);
  cen = max(range);

  % insert size chosen from Gaussian(s)
  isz = nan(nrows,1);
  pos=1;
  for i=1:size(distrib,1);
    nsub = round(distrib(i,1)*nrows);
    isz(pos:pos+nsub-1) = (distrib(i,2)+distrib(i,3)*randn(nsub,1));
    pos = pos+nsub;
  end

elseif ischar(distrib)

  % empirical distribution
  fprintf('Loading empirical insert-size distribution\n');
  distrib = load_isz_file(distrib);

end

if isstruct(distrib)

  T = distrib;
  require_fields(T,{'dat','mx'});

  % sum across lanes
  a = sum(T.dat,2);

  % if there's no data, return empty result
  if sum(a)==0, fprintf('No isz data!'); m=[]; return; end

  fprintf('Calculating expected stdev distribution\n');
  
  % choose cen (theoretical breakpoint location) = place where cumulative p > 0.99999
  ca = cumsum(a)/sum(a);
  cen = find(ca>0.99999,1)-1;
  if isempty(cen), cen = 2*T.mx; end

  % insert size chosen from empirical distribution
  a = round(a*nrows/sum(a));
  nrows = sum(a);
  isz = nan(nrows,1);
  pos=1;
  for i=0:T.mx
    isz(pos:pos+a(i+1)-1) = i;
    pos = pos + a(i+1);
  end

end

% choose starting point at random
x1 = cen*rand(nrows,1);
% add insert size
x2 = x1 + isz;

% filter reads:

% (old) keep all reads overlapping the center point (i.e. breakpoint)
%idx = find(x1<=cen & x2>=cen);

% (new) keep only those pairs with both reads aligned entirely on one side
idx = find(x1+P.readlen-P.maxpeek<=cen & x2-P.readlen+P.maxpeek>=cen);

% keep the filtered reads
x1 = x1(idx); x2 = x2(idx);
x1 = round(x1); x2 = round(x2);
nrows = length(idx);

% take permutations to calculate distribution of expected stdevs
%p = nan(length(setsizes),P.maxstdev);
m = zeros(length(setsizes),5);
for ss=1:length(setsizes)%,disp(setsizes(ss))
  q = x2(ceil(nrows*rand(setsizes(ss),P.nperms)));
  s = std(q);
  r = (max(q)+P.readlen-1)-min(q);
  m(ss,:) = [setsizes(ss) mean(s) std(s) mean(r) std(r)];
%  h = histc(s,1:P.maxstdev);
%  p(ss,:) = cumsum(h)/P.nperms;
end

return

figure(1),clf,hold on
plot(p')
xlim([0 200])
patch([30 30 100 100 30],[0 1 1 0 0],[1 1 0.8],'linestyle','none')
patch([40 40 90 90 40],[0 1 1 0 0],[1 1 0.6],'linestyle','none')
patch([50 50 80 80 50],[0 1 1 0 0],[1 1 0],'linestyle','none')
plot(p')
legend(num2str(setsizes'),'location','southeast')
xlabel('stdev of pos1');
ylabel('cumulative prob');
title('expected scatter of weirdpair sets');
hold off

keyboard

figure(2)
plot(log10(p'))
xlim([0 200])
legend(num2str(setsizes'),'location','southeast')
xlabel('stdev of pos1');
ylabel('log10 cumulative prob');
title('expected scatter of weirdpair sets');


