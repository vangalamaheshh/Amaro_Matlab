function plot_smoothed(CL21,si,cyto,fs,add_density,cols,lw)

r=[-1.5 1.5];
%lw=0.1;
%lw1=1;
%lw2=0.01;

if ~exist('cols','var')
  cols=[ 0 0 1; 1 0 0; 1 0 0];
end

if ~exist('add_density','var')
  add_density=1;
end

if ~exist('lw','var') || isempty(lw)
  lw.grid=0.1;
  lw.horizontal=1;
  lw.vertical=0.01;
end
  
if ~exist('fs','var') || isempty(fs)
  fs.title=5;
  fs.axes=5;
  fs.x_label=3;
  fs.y_label=3;
  fs.text=3;
end

if ~isstruct(fs)
  tmp.title=fs;
  tmp.axes=fs;
  tmp.x_label=fs;
  tmp.y_label=fs;
  tmp.text=fs;
  fs=tmp;
end

axis([0.5 size(CL21.dat,1)-0.5 r]);
axis on;
box on;
enlarge_axis([0 0.01]);
ax=axis;
chpos=find(diff(CL21.chrn));

CL21=add_cyto(CL21,cyto);

if (0)
  if mean(CL21.raw(:,1))>1
    disp('taking log of CL21.raw');
    CL21.raw(CL21.raw<0.1)=0.1;
    CL21.raw=log2(CL21.raw)-1;
  end
end

if ~isfield(CL21,'cbs_rl')
  CL21.cbs_rl{si}=runlength(CL21.dat(:,si),CL21.chrn);
  CL21.cbs_rl{si}=[ CL21.cbs_rl{si} CL21.chrn(CL21.cbs_rl{si}(:,1))];
end

rl=CL21.cbs_rl{si};
if (0)
  rl=rl(rl(:,4)~=23,:);
  if isfield(CL21,'medians')
    med=CL21.medians(si);
    rl(:,3)=rl(:,3)-med;
  else
    tmp=derunlength(rl);
    med=median(tmp);
    rl(:,3)=rl(:,3)-med;
  end
end

if add_density
  Nx=500;
  Ny=500;
  xstep=(ax(2)-ax(1))/(Nx-1);
  ystep=(ax(4)-ax(3))/(Ny-1);
  xv=(ax(1)-xstep):xstep:(ax(2)+xstep);
  yv=(ax(3)-ystep):ystep:(ax(4)+ystep);
  
  X=hist2d(1:size(CL21.raw,1),CL21.raw(:,si)-med,...
           xv,yv);
  X=X(2:(end-1),2:(end-1));
  %ph=plot(1:size(CL21.raw,1),CL21.raw(:,si),'.');
  %set(ph,'Color',[0.8 0.8 0.8]);
  %axis(ax);
  
  imagesc(xv(2:(end-1))+xstep/2,yv(2:(end-1))+ystep/2,X');
  set(gca,'YDir','normal');
  colormap(1-0.5*repmat((0:(1/63):1)',1,3));
  %keyboard
  %a1=gca;
  %a2=axes;
  %axis([ax -1 -0.5]);
  %imagesc(flipud(X'));
  %caxis([0 max(max(X(2:(end-1),2:(end-1))))]);
  %colormap(1-repmat((0:(1/63):1)',1,3));
end

armpos=find(diff(CL21.armn)>0);
for i=1:length(armpos)
  alh(i)=line([armpos(i) armpos(i)]+0.5,ax(3:4),'LineStyle','-','Color',[0.6 0.6 0.6],'LineWidth',lw.grid);  
end
for i=1:length(chpos)
  clh=line([chpos(i) chpos(i)]+0.5,ax(3:4),'LineStyle','-','Color',[0.3 0.3 0.3],'LineWidth',lw.grid);  
end
xlh=line(ax(1:2),[0 0],'LineStyle','-','Color',[0.3 0.3 0.3],'LineWidth',lw.grid);
if isfield(CL21,'sis')
  sample_name=CL21.sis(si).name;
else
  sample_name=CL21.sdesc(si);
end

if isfield(CL21,'used_normals')
  ttlh=title([sample_name ' (norm by ' num2str(length(CL21.used_normals{si})) ')'],'FontSize',fs.title);
else
  ttlh=title( sample_name,'FontSize',fs.title);
end
set(ttlh,'Interpreter','none');

set(gca,'XTick',[]);
set(gca,'FontSize',fs.axes);
chtxt=round(0.5*([ 0; chpos]+[chpos; size(CL21.dat,1)]));
th=text(chtxt,repmat(ax(3),length(chtxt),1)-0.04*(1.5+0.5*(-1).^(1:length(chtxt))')*(ax(4)-ax(3)),...
        CL21.chr(chtxt),'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',fs.x_label);


if (0)
  x=CL21.dat(:,si);
  plot(x,'k.'); hold on;
  y=clip_to_range(x,r);
  plot(find(x>r(2)),y(x>r(2),si),'r.');
  plot(find(x<r(1)),y(x<r(1),si),'b.');
end
if (0)
  rl=runlength(x);
  qq=CL21.cbs_rl{si};
  qq=qq(qq(:,4)~=23,:);
  range(qq(:,1:3)-rl(:,1:3))
end

%cols=[ 0 0 0; 1 0 0; 0 0 1];
rlh=[];1
for i=1:size(rl,1)
  if rl(i,3)>r(2)
    c=cols(2,:);
  elseif rl(i,3)<r(1)
    c=cols(3,:);
  else
    c=cols(1,:);
  end
  tmp1=line([rl(i,1)-0.5 rl(i,2)+0.5],clip_to_range(rl(i,3),r)*ones(1,2),...
              'LineWidth',lw.horizontal,'Color',c);
  if rl(i,1)~=rl(i,2)
    tmp1(2)=line([rl(i,1)-0.5 rl(i,1)+0.5],clip_to_range(rl(i,3),r)*ones(1,2),...
              'LineWidth',lw.horizontal,'Color',c);
    tmp1(3)=line([rl(i,2)-0.5 rl(i,2)+0.5],clip_to_range(rl(i,3),r)*ones(1,2),...
              'LineWidth',lw.horizontal,'Color',c);
  end
  if (i<size(rl,1)) && (rl(i,4)==rl(i+1,4))
    tmp2=line([rl(i,2)+0.5 rl(i,2)+0.5],clip_to_range([rl(i,3) rl(i+1,3)],r),...
              'LineWidth',lw.vertical,'Color',c); % [0.8 0.8 1]   
  elseif (i<size(rl,1))
    tmp2=line([rl(i,2)+0.5 rl(i,2)+0.5 rl(i,2)+0.5],clip_to_range([rl(i,3) 0 rl(i+1,3)],r),...
              'LineWidth',lw.vertical,'Color',c);
  else
    tmp2=[];
  end
  rlh=[rlh tmp1 tmp2];
end
%theory in practice
%2 stage


