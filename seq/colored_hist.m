function colored_hist(x,y,nbins);

if ~exist('nbins','var')
    [hy,hx]=hist(x);
else
    [hy,hx]=hist(x,nbins);
end    

for i=1:length(x)
    [mn,mni]=min(abs(hx-x(i)));
    xbin(i)=mni;
end
    

figure(1);
%subplot(2,1,1);
bar(hx,hy);
%subplot(2,1,2);
d=(hx(2)-hx(1))/2;
hgt=zeros(length(hx),1);
colormap('default')
for i=1:length(x)
    p = patch([hx(xbin(i))-d hx(xbin(i))+d hx(xbin(i))+d hx(xbin(i))-d],...
        [hgt(xbin(i)) hgt(xbin(i)) hgt(xbin(i))+1 hgt(xbin(i))+1],y(i,:));
%    set(gca,'CLim',[-2 2])
%    cdata = y(i, 1);
%    set(p,'FaceColor','flat',...
%          'FaceVertexCData',cdata,...
%          'CDataMapping','direct')
    hgt(xbin(i))=hgt(xbin(i))+1;
end
ylim([0 max(hy)+1]);
