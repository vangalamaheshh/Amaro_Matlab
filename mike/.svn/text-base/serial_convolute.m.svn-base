function S = serial_convolute(score,prob,n_iter,binsize)
%
% serial_convolute(score,prob,n_iter,binsize)
%
% returns a series of convolutions of the given histogram
%
% the histogram is enumerated as a matched set of scores and probabilities
% n_iter is the number of serial convolutions to perform
%
% default binsize = 0.05
%

if ~exist('binsize','var'), binsize = 0.05; end

maxscore = max(score) * n_iter;

% set up histogram

numbins = ceil(maxscore / binsize) + 1;   % bin1 = zero only
H = zeros(numbins,1);

% initial condition: all probability is in first bin (score=0, P=1)

H(1) = 1;

% sequential convolution

S = cell(n_iter,1);

prob(prob==0) = NaN;   % to avoid taking log of zero
offset = ceil(score/binsize);
offset = min(offset, numbins);

newH = zeros(numbins,1);

for i = 1:n_iter
  if ~mod(i,10), fprintf('%d/%d ',i,n_iter); end
  for s = 1:length(score)
    newH(:,s) = [repmat(0,offset(s),1) ; prob(s) * H(1:end-offset(s))];
  end
  H = sum(newH,2);
  S{i} = H;
end
fprintf('\n');


