function draw_primerplot(P,width)
% for analyzing Primer3 output
%
% Mike Lawrence 2010-04-06

if ~exist('width','var')
  seqlen = cellfun('length',P.SEQUENCE_TEMPLATE);
  width = max(seqlen);
end

Q = ones(slength(P),width);
for i=1:slength(P)
  if isnan(P.leftpos(i)) || isnan(P.rightpos(i))
    Q(i,:)=0;
  else
    Q(i,P.ampstart(i):P.ampend(i)) = 2;
    Q(i,P.leftpos(i):P.leftpos(i)+P.leftlen(i)-1) = 3;
    Q(i,P.rightpos(i)-P.rightlen(i)+1:P.rightpos(i)) = 3;
  end
end
imagesc(Q);
