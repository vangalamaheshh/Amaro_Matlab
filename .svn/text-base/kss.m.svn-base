function gscore=kss(geneset,dat)

dat=dna_norm(dat)./sqrt(size(dat,2)-1);
ng=size(dat,1);
for i=1:ng
  if mod(i,1000)==0
    disp(i)
  end
  dist=sum(dat-repmat(dat(i,:),ng,1).^2,2);
  [delta,deltai]=kol_smir(dist(geneset),dist(setdiff(1:size(dat,1), ...
                                                    geneset)));
  gscore(i)=delta;
end

