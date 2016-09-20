function [qvs,Y,pvs,ampdel]=plot_arm_histograms(fname,CL21,cyto,collapse_type,rng,stp,ts,set_ymax,horiz)

if ~exist('collapse_type','var')
  collapse_type='median';
end

if ~exist('ts','var')
  ts=[stp/2 stp/2];
end

if ~exist('horiz','var') || isempty(horiz)
  horiz=0;
end

CL21=add_cyto(CL21,cyto);
CL21.chrarmn=CL21.armn+CL21.chrn*2;

Y=collapse_D(CL21,'chrarmn',collapse_type); % 'mode'

% use all samples
samples_to_use=1:size(Y.dat,2);

armst={'p arm','q arm'};
figure(1); clf;
if horiz
  gr=make_subplotgrid([1.5 1*ones(1,23)],[1 1 0.7],...
                      1*ones(1,25),[ 0.2 2 1 0.2],0.001,0.3);
else
  gr=make_subplotgrid([ 1 4 4],[0.2 1*ones(1,23)],...
                      [ 1 0.5 0.5 1],2*ones(1,25),0.3,0.5);
end

u=unique(Y.chrn);
for i=1:length(u)
  if horiz
    subplotgrid(gr,3,u(i)+1);
    draw_names_box([ num2chromosome(u(i))],'vert',8,1,0,0.5,{'HorizontalAlignment','Center'});  % 'Chr'
  else
    subplotgrid(gr,u(i)+1,1);    
    draw_names_box([ num2chromosome(u(i))],'vert',8,1,1,0,{'HorizontalAlignment','Right'});  % 'Chr'
  end
end
if horiz
  for k=1:2
    subplotgrid(gr,k,1);
    draw_names_box(armst{k},'vert',8,1,0,1,{'HorizontalAlignment','Right'});
  end
else
  for k=1:2
    subplotgrid(gr,1,k+1);
    draw_names_box(armst{k},'vert',8,1,0.5,0);
  end
end

xvals=rng(1):stp:rng(2);
% ymax=max(sum(abs(Y.dat)<stp/2,2));
for i=1:size(Y.dat,1)
  if Y.n_collpase(i)<50  % skip small arms
    pvs(i)=NaN;
    continue;
  end
  h=hist(Y.dat(i,samples_to_use),xvals);
  left_side=find(xvals<(0-ts(2)));
  right_side=find(xvals>(0+ts(1)));
  left=sum(h(left_side));
  right=sum(h(right_side));

  H=zeros(size(h,2),3);
  H(left_side,1)=h(left_side);
  H(right_side,2)=h(right_side);
  H(:,3)=h'-sum(H(:,1:2),2);
  mx(i)=max([h(left_side) h(right_side)]);
end
ymax=max(mx);
if exist('set_ymax','var')
  ymax=set_ymax;
end

for i=1:size(Y.dat,1)
  if Y.n_collpase(i)<50  % skip small arms
    pvs(i)=NaN;
    continue;
  end
  if horiz
    subplotgrid(gr,Y.armn(i),Y.chrn(i)+1);
  else
    subplotgrid(gr,Y.chrn(i)+1,Y.armn(i)+1);    
  end
  h=hist(Y.dat(i,samples_to_use),xvals);
  left_side=find(xvals<(0-ts(2)));
  right_side=find(xvals>(0+ts(1)));
  left=sum(h(left_side));
  right=sum(h(right_side));

  H=zeros(size(h,2),3);
  H(left_side,1)=h(left_side);
  H(right_side,2)=h(right_side);
  H(:,3)=h'-sum(H(:,1:2),2);
  disp([ i max(h(left_side)) max(h(right_side)) ymax])
  bh=bar(rng(1):stp:rng(2),H,1,'stacked');
  set(bh,'LineStyle','none','LineWidth',eps,'ShowBaseLine','off');
  
  pv=cointest(left,left+right,0.5);
  pvs(i)=pv;
  lpv=(min(-log10(pv),6))/6;
  
  grc=0.7;
  if left>right
    c{1}=grc*ones(1,3)+[-grc*lpv -grc*lpv (1-grc)*lpv];
    c{2}=grc*ones(1,3);
    ampdel(i)=2;
  else
    c{2}=grc*ones(1,3)+[(1-grc)*lpv -grc*lpv -grc*lpv];
    c{1}=grc*ones(1,3);
    ampdel(i)=1;
  end
  c{3}=0.5*ones(1,3);
  for k=1:3
    set(bh(k),'FaceColor',c{k});
  end
%  ph=cut_bar(ymax,ymax*0.5,ymax*0.3,ymax*0.1,[1 1 1]);
  ax=axis;
  axis([ax(1:2) ax(3) ymax*1.7]);
 
  
%  ax=axis;
%  axis([rng(1)-stp/2 rng(2)+stp/2 0 ymax]);
  set(gca,'FontSize',4);
  set(gca,'YTick',[]);
  box off
  axis off
  th(i)=text(rng(2)-(rng(2)-rng(1))*0.35,ymax*1.5,['p=' num2str(roundsig(pv,2))],'FontSize',5);
  tick_sz=0.1;
%  line_color=[ 0 0 0];
  line_color=[ 0.5 0.5 0.5]; 
  line([rng(1)-stp/2 rng(2)+stp/2],[0 0],'Color',line_color);
  line([rng(1) rng(1)],[0 -ymax*tick_sz],'Color',line_color,'Clipping','off');
  line([0 0],[0 -ymax*tick_sz],'Color',line_color,'Clipping','off');
  line([rng(2) rng(2)],[0 -ymax*tick_sz],'Color',line_color,'Clipping','off');
  if Y.chrn(i)==1
    text(rng(1),-ymax*tick_sz,num2str(rng(1)),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',5);
    text(0,-ymax*tick_sz,'0','VerticalAlignment','top','HorizontalAlignment','center','FontSize',5);
    text(rng(2),-ymax*tick_sz,num2str(rng(2)),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',5);
  end
end

non_nan=find(~isnan(pvs));
pvs1=pvs(non_nan);
qvs1=calc_fdr_value(pvs1);
qvs=nan(size(pvs));
qvs(non_nan)=qvs1;

for i=1:length(non_nan)
  x=qvs1(i);
  x=roundsig(x,1);
  if x<0.01
    x=10^round(log10(x));
  end
  st=num2str(x);
  fi=find(st=='e');
  if ~isempty(fi)
    st=st(fi:end);
  end
  st=regexprep(st,'-0','-');
  if strmatch(st,'0.001')
    st='e-3';
  end
  if non_nan(i)<=2
    set(th(non_nan(i)),'String',['q=' st]);
  else
    set(th(non_nan(i)),'String',st);
  end    
end


if ~isempty(fname)
  print_D(fname,{{'fig'},{'pdf'},{'png','-r180'}}); % ,horiz);
end


