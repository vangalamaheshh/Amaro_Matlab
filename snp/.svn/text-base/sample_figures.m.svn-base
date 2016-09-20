function res=sample_figures(C,idx,cyto,fs,ext,start_i)

if ~exist('ext','var')
  ext=[];
end

if ~exist('start_i','var')
  start_i=1;
end
res=[];

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

if length(idx)>1
  for si=1:length(idx)
    res{si}=sample_figures(C,idx(si),cyto,fs,ext,start_i);
  end
  return
else
  i=idx;
end

if iscell(C)
  C2=C{2};
  C=C{1};
end


close all

%% smoothed data 
figure(1); clf;
set(gcf,'visible','off');
plot_smoothed(C,i,cyto,fs);
set(gcf,'visible','off');
ylabel('log2ratio','FontSize',3);
ttl=get(get(gca,'Title'),'String');
set(get(gca,'Title'),'String',[ 'Sample' sprintf('%03d',i+start_i-1) ': ' ttl]);
print_D(['Sample' sprintf('%03d',i-start_i+1) ext '.smoothed'],{{'fig'},{'png'}})

%% histogram analysis
figure(2); clf;
set(gcf,'visible','off');
rng=C.hist_rng;
midrng=(rng(1:(end-1))+rng(2:end))/2;
plot(midrng,C.spikes(:,i)/max(C.spikes(:,i))*max(C.hists(:,i)),'-k'); hold on
plot(midrng,C.hists(:,i)); hold on
set(gcf,'visible','off');
set(gca,'FontSize',fs.axes);
ax=axis;
axis([ midrng([min(find(C.hists(:,i)>1e-3*ax(4))) max(find(C.hists(:,i)>1e-3*ax(4)))]) ax(3:4)]);
enlarge_axis([0.1 0]);
ax=axis;

%% add lines for peak levels
for j=1:length(C.peaks{i})
  plot(C.peaks{i}(j),ax(4)*0.9,'r.','MarkerSize',3); 
end

for j=1:length(C2.peaks{i})
  plot(C2.peaks{i}(j),ax(4)*0.8,'g.','MarkerSize',3); 
end


if (0)
  notX=find(CL21.chrn~=23);
  r=CL21.raw(:,i)-CL21.medians(i);
  r=r(notX);
  [li,D,mv]=merge_levels(r,runlength(CL21.dat(notX,i)),0.5);% -log(1e-10)
  m=[];
  v=zeros(size(r,1),1);
  ax=axis;
  for j=1:length(li)
    m(j)=median(r(li{j}));
    v(li{j})=m(j);
    line([m(j) m(j)],ax(3:4),'Color','g');
  end
  [sm,smi]=sort(m);
  li=li(smi);
  %   2-(2.^(sm+1))
  %   2.^(sm+1)-2
  f1=find(sm>-1.3 & sm<-0.1);
  if ~isempty(f1)
    f1v=max(sm(f1));
    pure1=2-(2.^(f1v+1));
  else
    pure1=NaN;
  end
  f3=find(sm<0.9 & sm>0.1);
  if ~isempty(f3)
    f3v=min(sm(f3));
    pure3=2.^(f3v+1)-2;
  else
    pure3=NaN;
  end
  title(num2str([pure1 pure3]),'FontSize',fs);
end
set(gcf,'visible','off');
print_D(['Sample' sprintf('%03d',i-start_i+1) ext '.hist'],{{'fig'},{'png'}});


%% merged levels (on top of smoothed data)
figure(1);
set(gcf,'visible','off');
hold on
plot(C.dat(C.chrn~=23,i),'r-');
plot(C2.dat(C.chrn~=23,i),'g-');
set(gcf,'visible','off');
print_D(['Sample' sprintf('%03d',i-start_i+1) ext '.levels'],{{'fig'},{'png'}});  

%% pie charts
hv=histc(C2.level(:,i),1:length(C2.peaks{i}));
hv=as_row(hv)./sum(hv);
ueploid=C2.trans{i}.euploid;
if isempty(ueploid)
  up(i,1)=NaN;
  up(i,2)=NaN;
else
  up(i,1)=ueploid;
  up(i,2)=C2.peaks{i}(up(i,1));
end
close all
bw=0.01;
fse=4;
figure(1); clf
set(gcf,'visible','off');
bh=bar(C2.peaks{i},hv,0.1);
set(gcf,'visible','off');
ax=axis;
set(gcf,'visible','off');
pie_w=0.1*(ax(2)-ax(1));
tmp=get(gca,'DataAspectRatio');
pie_h=pie_w/tmp(1);
axis([ax(1)-pie_w/2 ax(2)+pie_w/2 -pie_h/2 ax(4)+pie_h/2]);
set(gcf,'visible','off');
enlarge_axis([ 0 0; 0 0.15]);
set(gcf,'visible','off');
ax=axis;
set(gcf,'visible','off');
tmp=get(gca,'DataAspectRatio');
pie_h=pie_w/tmp(1);
pp={};
light_color=[0.8 0.8 0.8];
delete(bh);
set(gcf,'visible','off');
line(ax(1:2),[0 0],[-eps -eps],'Color','k');
set(gcf,'visible','off');
if (1)
  for bi=1:length(C2.peaks{i})
    set(gcf,'visible','off');
    %    ph(bi)=patch([sm(bi)-bw/2 sm(bi)+bw/2 sm(bi)+bw/2 sm(bi)-bw/2 sm(bi)-bw/2],...
    %                 [0 0 ax(4) ax(4) 0],light_color);
    %    set(ph(bi),'EdgeColor','none');
    if ~isnan(up(i,1))
      text(up(i,2),hv(up(i,1))+pie_h/2+ax(4)*0.02,'Euploid','FontSize',fse,'HorizontalAlignment','center','VerticalAlignment','bottom');
      set(gcf,'visible','off');
    end
    sm=C2.peaks{i};
    phl(bi)=patch([sm(bi)-bw/2 sm(bi)+bw/2 sm(bi)+bw/2 sm(bi)-bw/2 sm(bi)-bw/2],...
                  [0 0 hv(bi) hv(bi) 0],[0 0 0]);
    set(gcf,'visible','off');
    set(phl(bi),'EdgeColor','none');
    set(gcf,'visible','off');
    a1=gca;
    set(gcf,'visible','off');
    
    if isfield(C2,'peaks_allele_balance')
      ab=C2.peaks_allele_balance{i};
    else
      ab=ones(1,length(C2.peaks{i}))
    end
    
    ab(ab==0)=eps;
    figure(2); clf;
    set(gcf,'visible','off');
    pp{bi}=pie([ab(bi),1-ab(bi)]);
    
    set(pp{bi}(1),'EdgeColor',[ 0 0 0]);
    set(gcf,'visible','off');
    set(pp{bi}(1),'FaceColor',[ 0 0 0]);
    set(gcf,'visible','off');
    if length(pp{bi}>2)
      set(pp{bi}(3:2:end),'FaceColor',light_color);
    end
    set(gcf,'visible','off');      
    
    figure(1);
    set(gcf,'visible','off');
    yv=hv(bi);
    put_in_3d(pp{bi}(1:2:end),...
              [sm(bi)-pie_w/2 sm(bi)+pie_w/2 sm(bi)+pie_w/2 sm(bi)-pie_w/2 sm(bi)-pie_w/2; ...
               yv-pie_h/2 yv-pie_h/2 yv+pie_h/2 yv+pie_h/2 yv-pie_h/2; zeros(1,5)]',0,1);
    set(gcf,'visible','off');      
  end
  close(2);
end
set(gcf,'visible','off');
print_D(['Sample' sprintf('%03d',i-start_i+1) ext '.barpie'],{{'fig'},{'png'}});

%% back to copy numbers


figure(3); clf;
set(gcf,'visible','off');
axis([-0.1 3.1 -0.1 3.1]);
colors='brgkm';
shapes='+xv';
r=2.^(C2.peaks{i}+1);
lh=[];
lg={};
trans=C2.trans{i};
cand=C2.trans{i}.CAND;
for cani=1:size(trans.P,1)
  for cei=1:size(trans.P,2)
    pval=squeeze(trans.PVAL(cani,cei,:));
    cn=squeeze(trans.CN(cani,cei,:));
    if nnz(pval<-0.5 | pval>1.5)==0 % should be 0 and 1
        lh(end+1)=plot(r,cn,[colors(cani) '-' shapes(cei)]); hold on;
        plot(2.^(C2.peaks{i}(cand(cani))+1),cn(cand(cani)),[colors(cani) 'o']); hold on;
        text(2.^(C2.peaks{i}(cand(cani))+1)+0.1,cn(cand(cani)),num2str(trans.P(cani,cei),'%3.2f'),'HorizontalAlignment','Left');
        lg{end+1}=[ num2str(pval)]; 
    end
  end
end
legend(lh,lg,'Location','Best');
ax=axis;
for level=ceil(ax(3)):floor(ax(4))
  line(ax(1:2),[level level],[-eps -eps],'Color',[0.8 0.8 0.8]);
end
set(gcf,'visible','off');
print_D(['Sample' sprintf('%03d',i) '.cn'],{{'fig'},{'png'}});


return



%% 
figure(3); clf;
set(gcf,'visible','off');
cloh=collapse_dat(L2t.dat(:,i),v,sm,'mean');
bar(sm,cloh,0.1);
print_D(['Sample' sprintf('%03d',i-start_i+1) '.loh'],{{'fig'},{'png'}});  
