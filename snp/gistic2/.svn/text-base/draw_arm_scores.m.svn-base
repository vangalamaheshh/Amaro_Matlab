function draw_arm_scores(D,arm_names,q,qv_cutoff,flip_direction,arm_name_side,bar_color,xlim_minqv)

if arm_name_side==-1 % left of figure
  disp_p=struct('x',struct('sizes',[0.1 0.03 0.15 0.5] ,'gaps',[1 1 1 0 1],'border',0.1),...
                'y',struct('sizes',[1],'gaps',[1 1],'border',0.2),...
                'items',...
                {{...
                    {1,1,'chrn','vert',6,-0.5},...
                    {1,2,'chrcyto'},...
                 }});
else % right of figure
  disp_p=struct('x',struct('sizes',[0.1 0.03 0.5 0.15] ,'gaps',[1 1 1 0 1],'border',0.1),...
                'y',struct('sizes',[1],'gaps',[1 1],'border',0.2),...
                'items',...
                {{...
                    {1,1,'chrn','vert',6,-0.5},...
                    {1,2,'chrcyto'},...
                 }});
end  

figure
[res,gr]=display_D(D,[],[],disp_p);

mid=zeros(length(arm_names),1);
start=mid;
stop=mid;
for i=1:length(arm_names)
  ci=str2num(arm_names{i}(1:(end-1)));
  if regexprep(arm_names{i}(end),' ','')=='p'
    ai=1;
  else
    ai=2;
  end
  disp([arm_names{i}  ' ' num2str(ci) ' ' num2str(ai)]);
  idx=find(D.chrn==ci & D.armn==ai);
  mid(i)=floor(mean(idx));
  start(i)=min(idx);
  stop(i)=max(idx);
end

if arm_name_side==-1
  subplotgrid(gr,1,4);
else
  subplotgrid(gr,1,3);
end

h=semilogy(mid,q,'o');
if exist('xlim_minqv','var')
    ylim([xlim_minqv 1]);
end
ax=axis;
delete(h);
set(gca,'YScale','log');
axis([1 length(D.chrn) ax(3:4)]);
set(gca,'XDir','reverse');
set(gca,'CameraUpVector',[-1 0 0]);
set(gca,'XTick',[],'tickdir','out');

% draw grey shading for even-numbered chromosomes and dotted lines for centromeres
for ci=1:max(D.chrn)
  inchr=find(D.chrn==ci);
  if ~isempty(inchr)
    mnc=min(inchr);
    mxc=max(inchr);
    if mod(ci,2)==0
      ph=patch([ mnc-0.5 mxc+0.5 mxc+0.5 mnc-0.5 mnc-0.5],...
               [ ax(3) ax(3) ax(4) ax(4) ax(3)],[ 0.9 0.9 0.9]);
      set(ph,'FaceAlpha',0.9,'EdgeColor','none');
    end
    cen=inchr(find(D.armn(inchr)==2,1));
    if ~isempty(cen)
      lh=line([cen-0.5 cen-0.5],[ax(3) ax(4)],'LineStyle',':');
      set(lh,'Color',[0.3 0.3 0.3]);
    end
  end
end
box on
% draw another box around the whole plot area to fix the gaps introduced by the grey shading
ax=axis;
line([ ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'Color',[0 0 0]);

tot=length(D.chrn);
prc=0.002;
for i=1:length(mid)
%  patch([mid(i)-tot*prc mid(i)+tot*prc mid(i)+tot*prc mid(i)-tot*prc],[1 1 q(i) q(i)],[1 0 0]);
  patch([start(i)+tot*prc stop(i)-tot*prc stop(i)-tot*prc start(i)+tot*prc],[1 1 q(i) q(i)],bar_color);
end
h=line([1 tot],[qv_cutoff qv_cutoff]);
set(h,'Color','g');
h=line([1 tot],[1 1]);
set(h,'Color','k');
if flip_direction
    camera_position=get(gca,'CameraPosition');
    camera_position(end)=-camera_position(end);
    set(gca,'CameraPosition',camera_position);
end

if arm_name_side==-1
  subplotgrid(gr,1,3);
else
  subplotgrid(gr,1,4);
end
sig_idx=find(q<qv_cutoff);
if arm_name_side==-1
  draw_spaced_text(arm_names(sig_idx),mid(sig_idx),1,length(D.chrn),0.1,0.3,0.32,0.5,0,[],...
                   'FontSize',9,'HorizontalAlignment','Right');
  set(gca,'XDir','reverse');
else
  draw_spaced_text(arm_names(sig_idx),mid(sig_idx),1,length(D.chrn),0.1,0.3,0.32,0.5,0,[],...
                   'FontSize',9);
end

