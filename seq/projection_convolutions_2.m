function p = projection_convolutions_2(Sdeg,Pdeg,score_obs,numbins,pmax,H,newH)
% Sdeg is (np,ncat+1,ncat+1)
% Pdeg is (np,ncat+1,ncat+1)
% score_obs is what p-value will be calculated against
% numbins is typically 1000
% pmax: once we are sure the pvalue is above this, we stop computing
% H and newH are allocated once in Matlab

if nargin~=7, error('wrong nargin'); end

np = size(Sdeg,1);
ncat = size(Sdeg,2)-1;
binsize = score_obs / numbins;
offset = min(numbins, round(Sdeg/binsize));
ncols = (ncat+1)*(ncat+2)/2;

if ~exist('H','var')
  H = zeros(numbins,1);
end
if ~exist('newH','var')
  newH = zeros(numbins,ncols);
end

H(:) = 0; H(1) = 1;  % initial condition: all probability is in first bin

% sequential convolution 
for p=1:np
  newH(:) = 0;
  col=1;
  for d1=0:ncat, for d2=0:d1
    o = offset(p,d1+1,d2+1);
    newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
    col=col+1;
  end,end
  H = sum(newH,2);
  p = max(0,1-sum(H));
  if (p>pmax) break; end
end



