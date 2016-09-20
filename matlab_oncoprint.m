%draft comut for showing mutation and deletion status. 

%generating toy data
A.TP53_del=[1 1 2 2 2 1 1 1 1 1 0 0 0];
A.TP53_mut=[1 1 0 0 0 1 1 1 1 1 0 1 0];
A.TP53_methyl=[0 0 0 0 0 0 0 0 0 0 0 1 0];

A.MYC_amp=[1 0 0 0 0 0 0 0 1 1 1 2 3];

%number of variants
NV=2;

%number of samples
N=length(A.TP53_del);


%generating colors for legend
hp=plot(1,NV,'s',1,NV,'s',1,NV,'s');
set(hp(1),'markerfacecolor',0.5*[0 0 .1],'color',0.5*[0 0 .1]);
set(hp(2),'markerfacecolor',0.4*[1 1 1],'markeredgecolor',0.5*[0 0 .1]);
set(hp(3),'markerfacecolor',0.4*[1 1 1],'markeredgecolor',0.5*[0 0 .1]);



%setting axis
axis([0 N+0.5 0.5 NV+0.5]); set(gca,'visible','off');
set(gca,'box','on')


dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;


for s=1:N
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
if A.TP53_del(s)>0
patch(x1,y1,(1/A.TP53_del(s))*[65/255 105/255 225/255],'edgecolor',[.9 .9 .9])
else
    if A.TP53_methyl(s)==1
       patch(x1,y1,[.9 .9 .9],'edgecolor',[65/255 105/255 225/255])
    else
patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])
    end
end
end

for s=1:N
x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
    j=NV;
    y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
    
    if A.TP53_mut(s)==1
    patch(x1,y1,[255/255 215/255 0/255],'edgecolor',[255/255 215/255 0/255])
    end

end

for s=1:N
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-1;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
if A.MYC_amp(s)>0
patch(x1,y1,(1/A.MYC_amp(s))*[255/255 105/255 65/255],'edgecolor',[.9 .9 .9])
else
    
patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])
end
end


lab{1}='TP53';
lab{2}='ONC';
text(xlim(1)-1,NV,lab{1},'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',16)
text(xlim(1)-1,NV-1,lab{2},'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',16)

set(gca,'ytick',1:NV,'yticklabel',fliplr(lab))
set(gca,'xticklabel',[])
box off




