function my_imagesc(dat,sv)

mn=min(dat(:));
mx=max(dat(:));

m=1+sv+floor((dat-mn)/(mx-mn)*(64-sv));
image(m);
caxis([mn mx]);
h=colorbar;
set(get(h,'Children'),'CData',((sv+1):64)');
keyboard
