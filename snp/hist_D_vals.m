function hist_D_vals(T,rng,supid)

hist(T.dat(~isnan(T.dat(:,1)),1),rng);
ax=axis; axis([min(rng) max(rng) ax(3:4)]);
idx=min(find(~isnan(T.dat(:,1))));
armst='pq';
st=[T.chr{idx} armst(T.armn(idx))];
cols='rgky';
for i=1:length(supid)
  line(ones(1,2)*T.supdat(supid(i),1),ax(3:4),'Color',cols(mod(i-1,length(cols))+1));
  st=[st ' ' num2str(T.supdat(supid(i),1))];
end
title(st);

