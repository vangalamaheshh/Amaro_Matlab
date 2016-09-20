function H = batch_convolute(Y)
%
% batch_convolute(Y)
%
% returns the convolutions of the given histograms
%
% each cell of Y is an expanded-format histogram;
% all are assumed to have the same binsize.
%

ny = length(Y);

numbins = 1;
for i=1:ny
  if size(Y{i},1)==1, Y{i}=Y{i}'; end   % make sure all are row vectors
  numbins = numbins + length(Y{i})-1;
end

% set up histogram

H = zeros(numbins,1);

% initial condition: all probability is in first bin (score=0, P=1)

H(1) = 1;

% method 2


for i=1:ny
%i
  newH = zeros(numbins,1);
  idxh = find(H>0);
  for hi=1:length(idxh)
%fprintf('%d/%d ',hi,length(idxh));
    hbin = idxh(hi);
    add = Y{i}*H(hbin);
    range = hbin+[1:length(add)]-1;
    newH(range) = newH(range) + add;
  end
  H = newH;
end



if 0

% method 1

for i=1:ny
  newH = zeros(numbins,1);
  idxh = find(H>0);
  idxy = find(Y{i}>0);
  for hi=1:length(idxh)
    for yi=1:length(idxy)
      hbin = idxh(hi);
      ybin = idxy(yi);
      newbin = (hbin-1 + ybin-1) + 1;
      newprob = H(hbin) * Y{i}(ybin);
      newH(newbin) = newH(newbin) + newprob;
    end
  end
  H = newH;
end

end
