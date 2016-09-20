function classify_dRanger_input(X)

NN=1; CHR1=2; STR1=3; ST1=4; END1=5;
CHR2=6; STR2=7; ST2=8; END2=9;
QUAL1=10; QUAL2=11; RGRP=12;

inter = X(:,CHR1)~=X(:,CHR2);
intra = ~inter;
strand_normal = (X(intra,STR1)==0 & X(intra,STR2)==1 & X(intra,ST1)<X(intra,ST2));
span = X(intra,ST2)-X(intra,ST1);
upper = [inf  10e6 1e6   100e3 10e3 5e3 2e3 1e3 600];
lower = [10e6 1e6  100e3 10e3  5e3  2e3 1e3 600 -inf];

fprintf('Inter\nIntra\n');
fprintf('\nIntra (strand-normal)\n');
for i=1:length(upper), fprintf('%d-%d\n',lower(i),upper(i)); end
fprintf('\nIntra (strand-weird)\n');
for i=1:length(upper), fprintf('%d-%d\n',lower(i),upper(i)); end
fprintf('\n');
fprintf('%0d\n%0d\n',sum(inter),sum(intra));
fprintf('\n\n');
for i=1:length(upper)
  fprintf('%0d\n',sum(span(strand_normal)<upper(i) & span(strand_normal)>=lower(i)));
end
fprintf('\n\n');
for i=1:length(upper)
  fprintf('%0d\n',sum(span(~strand_normal)<upper(i) & span(~strand_normal)>=lower(i)));
end


