function P = semiexactP(Pdistrib,Sobs,binsize)
%%
%% semiexactP(Pdistrib,Sobs,binsize)
%%
%% inputs:
%%    Pdistrib matrix (i rows, N+1 columns)
%%        probability distribution for each item
%%        #rows = number of distributions to be convoluted
%%        #cols = number of states in each distribution
%%        rows = items (distributions)
%%        cols = probabilities across possible states
%%   Sobs vector (i rows)
%%        observed state for each item (range = 0 to N)
%%   binsize
%%        # of bins per log10
%%
%% procedure:
%%   1. computes Pobs, the point probability of the given state
%%   2. convolutes distributions, discarding probabilities less likely than Pobs
%%
%% output:
%%   P, the probability of getting an outcome at least as extreme as the observed one.
%%
%% Mike Lawrence 2008-04-25
%%


% check input parameters

[ni ns] = size(Pdistrib);

if length(Sobs) ~= ni
   error('Sobs and Pdistrib must have same number of rows.');
end

last_state = ns-1;

if any(Sobs<0 | Sobs>last_state | Sobs~=round(Sobs))
   error('Sobs must be an integer between 0 and N, where N+1 is # of columns in Pdistrib');
end

% compute Pobs, the point probability of the observed state

Pobs = 1;
for i=1:ni
    Pobs = Pobs * Pdistrib(i,Sobs(i)+1);
end
if Pobs==1, P=1; return; end
if Pobs==0, Pobs=10^-20; end
score_obs = -log10(Pobs);

% set up histogram

numbins = ceil(score_obs / binsize) + 1;   % bin1 = zero only
H = zeros(numbins,1);

% initial condition: all probability is in first bin (score=0, P=1)

H(1) = 1;

% sequential convolution

Pdistrib(Pdistrib==0) = NaN;   % to avoid taking log of zero
offset = ceil(-log10(Pdistrib)/binsize);
offset = min(offset, numbins);

newH = zeros(numbins,ns);

for i = 1:ni
  for s=0:last_state
    o = offset(i,s+1);
    newH(o+1:end,s+1) = Pdistrib(i,s+1) * H(1:end-o); 
%    newH(:,s+1) = [repmat(0,offset(i,s+1),1) ; Pdistrib(i,s+1) * H(1:end-offset(i,s+1))];
  end
  H = sum(newH,2);
  newH(:) = 0;
end

% return result

P = 1 - sum(H);
