function r = score_methyl(ALL)

r = zeros(ALL.ng,1);
for i=1:ALL.ng
  if ~ALL.gene.have_methyl(i), continue; end
  x = 1-ALL.methyl(i,:);
  y = ALL.expr(i,:);
  good = find(~isnan(x) & ~isnan(y));
  x=x(good);
  y=y(good);
  tmp = corrcoef(x,y);
  r(i) = tmp(1,2);
  if isnan(r(i)), r(i)=0; end
  if abs(r(i))>0.5 && 0
    scatter(x1,y1);
    title([ALL.gene.name{i} ' r=' num2str(r(i))]);
    keyboard
  end
 end
