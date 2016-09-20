function [chr pos] = convert_wins(wins,windowsize)
if ~exist('windowsize','var'), windowsize=100;end

clen = load_chrlen;
nwins = round(1.1*clen/windowsize);

cs = cumsum(nwins);

chr = nan(length(wins),1);
pos = chr;

for i=1:length(wins)
  chr(i) = find(wins(i)<cs,1);
  if chr(i)==1, subtr=0; else subtr=cs(chr(i)-1); end
  pos(i) = (wins(i)-subtr)*windowsize;
end



