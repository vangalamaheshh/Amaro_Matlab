function [pmax pmin] = projection_1d_convolutions(Sdeg,Pdeg,score_obs,numbins,H,newH)
% Sdeg is (np,ncat+1), must be integers
% Pdeg is (np,ncat+1)
% score_obs a sum of integer values from Sdeg
% numbins is typically 2*score_obs to leave a margin for finding pmin
% H and newH are allocated once in Matlab

np = size(Sdeg,1);
ncat = size(Sdeg,2)-1;
ncols = ncat+1;

% binsize is fixed at 1

if ~exist('numbins','var'), numbins = max(10,score_obs*2); end
if ~exist('H','var'), H = zeros(numbins,1); end
if ~exist('newH','var'), newH = zeros(numbins,ncols); end

H(:) = 0; H(1) = 1;  % initial condition: all probability is in first bin

% sequential convolution 
for p=1:np
  newH(:) = 0;
  col=1;
  for d=0:ncat
    o = Sdeg(p,d+1);
    newH(o+1:end,col) = Pdeg(p,d+1) .* H(1:end-o);
    col=col+1;
  end
  H = sum(newH,2);
end

%sparse(H)

% p-value bounds
pmax=1;
pmin=1;
for i=1:length(H)
  if H(i)>0
    pmin=pmin-H(i);
    if i>score_obs, break; end
    pmax=pmax-H(i);
  end
end

if (pmax<0), pmax=0; end
if (pmin<0), pmin=0; end
