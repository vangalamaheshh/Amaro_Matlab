function draw_CNQC_plot(Lt,Ln,R,T,N,Z,individual_name)
% Mike Lawrence 2009-10-28

have_tumor_seg = ~all(isnan(Z(:,1)));
have_normal_seg = ~all(isnan(Z(:,2)));

I = []; c = []; n = [];
w = size(T,2) + size(N,2);spw = ceil(w/50);segw = ceil(w/8);sp = zeros(size(T,1),spw);

if have_tumor_seg
  I = [I repmat(Z(:,1),1,segw) sp];
  c = [c ones(1,segw) zeros(1,spw)];
  n = [n nan(1,segw+spw)]; 
end

% tumor lanes
I = [I T(:,Lt.enoughreads)];
c = [c 2*ones(1,sum(Lt.enoughreads))];
n = [n find(Lt.enoughreads)'];

if have_normal_seg
  I = [I sp repmat(Z(:,2),1,segw)];
  c = [c zeros(1,spw) 3*ones(1,segw)];
  n = [n nan(1,spw+segw)];
end

% normal lanes
I = [I sp N(:,Ln.enoughreads)];
c = [c zeros(1,spw) 4*ones(1,sum(Ln.enoughreads))];
n = [n nan(1,spw) find(Ln.enoughreads)'];

I(isnan(I))=1;

% restrict width of divider bars to 1/100

mincols = 100;
ncols = size(I,2);
if any(c==0) && ncols<mincols
  xf = ceil((mincols-ncols)/ncols);
  ord = reshape(repmat(1:ncols,xf,1),1,xf*ncols);
  del = find(c(ord(1:end-1))==0 & c(ord(2:end))==0);
  ord(del) = [];
  c = c(ord);
  I = I(:,ord);
  n = n(:,ord);
end

% display image

imagesc(I,[0 2]);
colormap(bwr);
colorbar;
ylabels_by_group(R.chr);
set(gca,'xtick',[]);title([individual_name ' CopyNumberQC'],'fontsize',20,'interpreter','none');
ty = 1.015*size(I,1);
txtparams = {'fontsize',18,'horizontalalign','center','verticalalign','top'};
if have_tumor_seg, text(mean(find(c==1)),ty,'S_T',txtparams{:}); end
text(mean(find(c==2)),ty,'T lanes',txtparams{:});
if have_normal_seg, text(mean(find(c==3)),ty,'S_N',txtparams{:}); end
text(mean(find(c==4)),ty,'N lanes',txtparams{:});

% make divider bars black instead of blue

for i=1:length(c)
  if c(i)==0
    rectangle('position',[i-0.5 -size(I,1)*.01 1 size(I,1)*1.1],...
      'facecolor',[0 0 0],'edgecolor',[0 0 0],'clipping','off')
  end
end

% draw red/pink asterisks under problem lanes, and blue asterisks under blacklisted lanes

yy = 1.03*size(I,1);
txtparams = {'fontsize',20,'horizontalalign','center','verticalalign','middle'};
% tumor mixups
mixups = grep('MIXUP|CONTAM',Lt.judgement,1);
for i=1:length(mixups)
  x = mean(find(c==2 & n==mixups(i)));
  text(x,yy,'*','color',[1 0 0],txtparams{:});
end  
% normal mixups
mixups = grep('MIXUP|CONTAM',Ln.judgement,1);
for i=1:length(mixups)
  x = mean(find(c==4 & n==mixups(i)));
  text(x,yy,'*','color',[1 0.53 0.53],txtparams{:});
end
% tumor blacklisted lanes
bl = find(Lt.is_blacklisted);
for i=1:length(bl)
  x = mean(find(c==2 & n==bl(i)));
  text(x,yy,'*','color',[0 0 1],txtparams{:});
end
% normal blacklisted lanes
bl = find(Ln.is_blacklisted);
for i=1:length(bl)
  x = mean(find(c==4 & n==bl(i)));
  text(x,yy,'*','color',[0.5 0.5 1],txtparams{:});
end
