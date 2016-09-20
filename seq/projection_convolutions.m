function p = projection_convolutions(Sdeg,Pdeg,score_obs,numbins,H,newH)
% Sdeg is (np,ncat+1,ncat+1)
% Pdeg is (np,ncat+1,ncat+1)
% numbins is typically 1000
% H and newH are allocated once in Matlab

np = size(Sdeg,1);
ncat = size(Sdeg,2)-1;
binsize = score_obs / numbins;
ncols = (ncat+1)*(ncat+2)/2;
offset = min(numbins, round(Sdeg/binsize));

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
%    fprintf('p=%d d1=%d d2=%d idx=%d o=%d\n',p,d1,d2,nan,o);
    newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
    col=col+1;
  end,end
  H = sum(newH,2);
end

% p-value
p = max(0,1-sum(H));



