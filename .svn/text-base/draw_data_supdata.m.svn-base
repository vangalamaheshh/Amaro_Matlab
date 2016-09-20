function draw_data_supdata(D,gidx,sidx,supidx,nanval)

if isstruct(D)
  dat=D.dat;
else
  dat=D;
end

if isempty(gidx)
  gidx=1:size(dat,1);
end

if isempty(sidx)
  sidx=1:size(dat,2);
end

dat=dat(gidx,sidx);
dat(isnan(dat))=nanval;
subplot(7,1,1:3);
imagesc(dat);
bluepink;
colorbar;
subplot(7,1,4:6);
imagesc_trim(dna_norm(dat));
bluepink;
colorbar;
subplot(7,1,7);
xx=D.supdat(supidx,sidx);
xx(isnan(xx))=nanval;
imagesc(dna_norm(xx));
bluepink;
colorbar;
