function draw_ttype_bars_figure(X,outname,P)

if ~exist('P','var'), P=[]; end

if ~grepmi('\.eps$',outname), error('outname must end with .eps'); end

% split long names
X.longname_split = X.longname;
X.longname_split{strcmp('Chronic lymphocytic leukemia',X.longname)} = {'Chronic lymphocytic','leukemia'};
X.longname_split{strcmp('Diffuse large B-cell lymphoma',X.longname)} = {'Diffuse large','B-cell lymphoma'};
X.longname_split{strcmp('Esophageal adenocarcinoma',X.longname)} = {'Esophageal','adenocarcinoma'};
X.longname_split{strcmp('Glioblastoma multiforme',X.longname)} = {'Glioblastoma','multiforme'};
X.longname_split{strcmp('Lung squamous cell carcinoma',X.longname)} = {'Lung squamous','cell carcinoma'};
X.longname_split{strcmp('Acute myeloid leukemia',X.longname)} = {'Acute myeloid','leukemia'};

% frequency symbol definitions
F=[];
F.cutoff = [1 2 3 5 10 inf]';
F.color = {[0.4 0 0.3],[0 0.3 0.9],[0 0.8 1],[0 0.7 0],[1 0.5 0],[1 0 0]}';
F.ss = [7 12 14 16 18 20]'/4000;   % symbol size
F.ssq = [6 3 2 1 0 -1]'/4000;      % symbol x-position fine adjustment

% tier colors
cols = [0 0.3 1;0 0.3 1;0 0.3 1;0 0 0;1 0 0];

% DRAW FIGURE

ti_to_show = grepvi('pancan|union|best_ttype|leaveout|\_|COAD|READ',X.name,1); % show only 21 tumor types
ntt = length(ti_to_show);
[tmp ord] = sort(X.longname(ti_to_show)); ti_to_show = ti_to_show(ord); % sort by longname

close all,figure(1),clf,set(gca,'visible','off','position',[0 0 1 1]);
set(gcf,'color',[1 1 1],'position',[400 -20 630 840]);hold on;ylim([0 1]);xlim([0 1]);
set(gcf,'Units','Inches'); pos = get(gcf,'Position'); set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ytop = 0.95; ybtm = 0.01; xname=0.20; xleft = 0.22; xright = 0.95;

for tii=1:ntt;
  ti=ti_to_show(tii);
  y = ytop - ((tii-1.3)/(ntt-1)*(ytop-ybtm));
  line([xleft xright],[y y],'color',[0 0 0]); npat = ['\fontsize{6}n=' num2str(X.npat(ti))]; ntxt = X.longname_split{ti};
  if ~iscell(ntxt), ntxt={ntxt}; end;
  ntxt = [ntxt npat]; text(xname,y,ntxt,'color',[0 0 0],'fontsize',10,'horizontalalignment','right','verticalalignment','bottom');

  % draw yellow and orange regions
  xpad = 0.002;
  numsig=sum(X.G{ti}.q<=0.1); gp=X.G{ti}.pmax(numsig+[0 1]); glp=min(15,-log10(gp)); glp=mean(glp); xthresh=xright-glp*(xright-xleft)/15;
  line([xleft-xpad xthresh],[y y],'linewidth',3,'color',[1 1 0]);
  numnear=20; nearidx = numsig+[1:numnear]; nearidx_featured = nearidx(ismember(X.G{ti}.tier(nearidx),[1 2 3 5])); gp = X.G{ti}.pmax(nearidx);
  glp = min(15,-log10(gp)); x = xright - glp*(xright-xleft)/15; line([xthresh max(x)+xpad],[y y],'linewidth',3,'color',[1 0.7 0]);

  % show all significant genes
  gidx = 1:numsig; gname = X.G{ti}.gene(gidx); gp = X.G{ti}.pmax(gidx); glp = min(15,-log10(gp)); x = xright - glp*(xright-xleft)/15;
  xadj_dot = find_text_pos(x,0.0005,0.0005,xleft,xright); xadj_txt = find_text_pos(x,0.005,0.005,xleft,xright);

  % y-offset definitions
  texty=0.011; doty=0.005; ly1=0.005; ly2=0.004; ly3=0.002; ly4=-0.002;

  % gene names
  ccc = X.G{ti}.tier(gidx);ccc(ccc==0)=4;
  h = textextra(xadj_txt,repmat(y+texty,length(xadj_txt),1),gname,'fontsize',5,'color',nansub(cols,ccc),'rotation',90,'horizontalalignment','left');

  % frequency symbols
  for i=1:numsig, pct=100*X.G{ti}.npat(i)/X.npat(ti); k=find(pct<=F.cutoff,1);
    fdot(xadj_txt(i)-0.75*F.ss(k), y+doty, k);
  end

  % draw lines connecting genes to dots
  p = {'color',[0.2 0.2 0.2],'linewidth',0.2};
  for i=1:length(x)
    line([xadj_txt(i) xadj_txt(i)],y+[ly1 ly2],p{:});
    line([xadj_txt(i) xadj_dot(i)],y+[ly2 ly3],p{:});
    line([xadj_dot(i) xadj_dot(i)],y+[ly3 ly4],p{:});
  end

  % over-the-horizon "next 20 genes"
  nnear = length(nearidx_featured);
  fmod = 4; if ismember(X.name{ti},{'CLL','DLBCL','MED'}), fmod=5; end  % how many names to list per line
  if ismember(X.name{ti},{'CRC'}), fmod=3; end

  fx0 = min(0.91,max(xthresh,max(xadj_txt)+0.02));
  fyspace = 0.009; fxspace = 0.0017;
  fy0 =  y + 0.008 + fyspace*ceil(nnear/fmod);
  if strcmp(X.name{ti},'UCEC'), fx=xthresh; fy=y-0.007; else fx=fx0; fy=fy0; end
  params = {'fontsize',5,'verticalalignment','middle'};

  ftxt =  ['next ' num2str(numnear) ' genes']; if length(nearidx_featured)>0, ftxt = [ftxt ' include']; end
  h = text(fx,fy,ftxt,params{:}); e = get(h,'extent');
  if strcmp(X.name{ti},'UCEC'), fx = fx + e(3) + fxspace; else fy=fy-fyspace; end

  for ni=1:nnear, i=nearidx_featured(ni);
    gname = X.G{ti}.gene{i}; tier = X.G{ti}.tier(i); if tier==0, tier=4; end
    pct=100*X.G{ti}.npat(i)/X.npat(ti); k=find(pct<=F.cutoff,1);
    fdot(fx, fy-F.ss(k)*0.55, k);
    fx=fx+F.ss(k)*1.5+0.001; % space between dot and name
    h = text(fx,fy,gname,params{:},'color',nansub(cols,tier)); e = get(h,'extent');
    if mod(ni,fmod)==0, fx=fx0; fy=fy-fyspace; else fx = fx + e(3) + fxspace; end
  end

end  % next tumor type
hold off

% p-value axis
for glp=15:-1:0
  py=0.007;x=xright-glp*(xright-xleft)/15;txt=['10^{-' num2str(glp) '}'];if glp==0,txt='1';elseif glp==1,txt='0.1';elseif glp==2,txt='0.01';end
  text(x,py,txt,'color',[0 0 0],'fontsize',6,'horizontalalignment','center'); line([x x],py+[0.007 0.012],'color',[0 0 0],'linewidth',0.5);
end, text(xname,py,'gene p-value','color',[0 0 0],'fontsize',6,'horizontalalignment','right');

% frequency legend
lx = 0.03; ly=0.07; lsy=0.012; ct=num2cellstr(F.cutoff);
for k=1:length(F.cutoff), y=ly+(k-1)*lsy;
  if k==1, ftxt = ['<' ct{k} '%']; elseif k==length(F.cutoff), ftxt = ['>' ct{k-1} '%']; else ftxt = [ct{k-1} '-' ct{k} '%']; end
  fdot(lx-0.75*F.ss(k), y-0.003, k);
  text(lx+0.015,y,ftxt,'fontsize',6);
end
text(lx,ly+k*lsy,'frequency','fontsize',6);
rectangle('position',[lx-0.01 ly-0.01 lx+0.05 ly+0.022]);

% print figure
print(outname,'-depsc');



   function fdot(x,y,k)
     rectangle('position',[x y F.ss(k)*1.5 F.ss(k)*1.1],'curvature',[1 1],'facecolor',F.color{k},'linestyle','none');
   end


end














