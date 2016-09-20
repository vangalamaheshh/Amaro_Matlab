function plot_scatter_regs(CL21,regs,ts,cyto,ext,r,m,rgg,use_regs_for_resid_qv,supid)
%% add support for supid and supmark
if ~exist('ext','var')
  ext=[];
end

if ~exist('r','var')
  r=[-1 1];
end

if ~exist('n','var')
  m=[0.1 0.1];
end

if ~exist('use_regs_for_resid_qv','var')
  use_regs_for_resid_qv=0;
end


y1=CL21.dat(:);
pct=[5 10 25 50 75 90 95];
pct_y1=prctile(y1,pct);
pct_z1=clip_to_range(pct_y1,r,m);

kstr={'amp','del'};

loc_chrn=cat(1,rgg.chrn);
loc_st=cat(1,rgg.start);
loc_en=cat(1,rgg.end);


for k=1:2
  sfname=['scatter_' kstr{k} ext '.ps'];
  if exist(sfname,'file')
    delete(sfname);
  end
  figure(1); clf;
  for i=1:length(regs{k})
    subplot(1,5,mod((i-1),5)+1);
    wpk(1)=regs{k}(i).peak_wide_st;
    wpk(2)=regs{k}(i).peak_wide_en;
    wpk=wpk(1):wpk(2);
    x=i+((1:size(CL21.dat,2))/size(CL21.dat,2)-0.5)*0.8;
    y=CL21.dat(regs{k}(i).peak,:);
    z=clip_to_range(y,r,m);

    if exist('supid','var')
      shape='^ov';
      for j=1:max(CL21.supdat(supid,:))
        idx=find(CL21.supdat(supid,:)==j);
        ph=plot(x(idx),z(idx),shape(j),'MarkerSize',3,'MarkerFaceColor','b'); hold on
      end
    else
      ph=plot(x,z,'o','MarkerSize',3,'MarkerFaceColor','b'); hold on      
    end
    axis([i-0.5 i+0.5 r(1)-m(1)*1.1 r(2)+m(2)*1.1]);
    ax=axis;
    yt=get(gca,'YTick');
    ytl=get(gca,'YTickLabel');
    inside=find(yt>=r(1) & yt<=r(2));
    yt=yt(inside);
    set(gca,'YTick',yt,'YTickLabel',deblank(cellstr(ytl(inside,:))));
    if mod(i-1,5)~=0
      set(gca,'YTickLabel',[]);
    end
    set(gca,'TickLength',[0 0 ]);
    for ti=1:length(yt)
     set(gca,'TickLength',[0 0 ]);
     line([ax(1) ax(1)+(ax(2)-ax(1))*0.05],[yt(ti) yt(ti)],'Color','k');
    end     
     
    if k==1
      above=find(y>ts(1));
    else
      above=find(y<-ts(2));
    end
    if exist('supid','var')
      shape='^ov';
      for j=1:max(CL21.supdat(supid,:))
        idx=find(CL21.supdat(supid,:)==j);
        plot(x(intersect(above,idx)),z(intersect(above,idx)),...
             [ 'r' shape(j)],'MarkerSize',3,'MarkerFaceColor','r'); hold on
        nabove_sup(j)=length(intersect(above,idx));
      end
    else
      plot(x(above),z(above),...
           'ro','MarkerSize',3,'MarkerFaceColor','r'); hold on      
    end
    
    nabove=length(above);
    plot(ax(1:2),ts(1)*ones(1,2),'k-');
 %   plot(ax(1:2),(log2(2.4)-1)*ones(1,2),'k--');
    plot(ax(1:2),-ts(2)*ones(1,2),'k-');
 %   plot(ax(1:2),(log2(1.6)-1)*ones(1,2),'k--');
    
    
    [ol,rng]=outliers(y,'iqr');
    plot(ax(1:2),[rng(1) rng(1)],'k:');
    plot(ax(1:2),[rng(2) rng(2)],'k:');
    
    
    %plot(x(ol),z(ol),'ro');
    
    [st,chr,bp]=genomic_location(CL21,{[regs{k}(i).peak_st regs{k}(i).peak_en]},cyto,1,1);
    [Wst,Wchr,Wbp]=genomic_location(CL21,{[regs{k}(i).peak_wide_st regs{k}(i).peak_wide_en]},cyto,1,1);
    gidx=find(loc_chrn==chr & abs((loc_st+loc_en)*0.5-mean(bp{1}(2),bp{1}(1)))<2e6);
    
    set(gca,'XTick',i);
    pct=[5 10 25 50 75 90 95];
    pct_y=prctile(y,pct);
    pct_z=clip_to_range(pct_y,r,m);
    ph(1)=patch(ax([ 1 2 2 1 1]),[r(1) r(1) r(1)-m(1)*1.1 r(1)-m(1)*1.1 r(1)],-ones(1,5),[0.9 0.9 0.9],'LineStyle','none');
    ph(2)=patch(ax([ 1 2 2 1 1]),[r(2) r(2) r(2)+m(2)*1.1 r(2)+m(2)*1.1 r(2)],-ones(1,5),[0.9 0.9 0.9],'LineStyle','none');
    line(ax([1 2 2 1 1]),ax([ 3 3 4 4 3]),'Color','k');
    box on
    
    for pi=1:length(pct_z)
      line([ax(2)-(ax(2)-ax(1))*0.05 ax(2)],[pct_z(pi) pct_z(pi)],'Color','k');
%      text(ax(2),pct_z(pi),num2str(pct(pi)),'FontSize',5);
    end
    for pi=1:length(pct_z1)
      line([ax(2)-(ax(2)-ax(1))*0.05 ax(2)],[pct_z1(pi) pct_z1(pi)],'Color','r');
      text(ax(2),pct_z1(pi),num2str(pct(pi)),'FontSize',5);
    end
    
    plot(ax(1:2),clip_to_range([pct_y1(3)-(pct_y1(5)-pct_y1(3)) pct_y1(3)-(pct_y1(5)-pct_y1(3))],r,m),'r-');
    plot(ax(1:2),clip_to_range([pct_y1(5)+(pct_y1(5)-pct_y1(3)) pct_y1(5)+(pct_y1(5)-pct_y1(3))],r,m),'r-');
    
    if use_regs_for_resid_qv && abs(log2(regs{k}(i).qv)-log2(regs{k}(i).resid_qv))<0.0001 % nearly identical
      resid_qv_st=[' resid. qv:' num2str(regs{k}(i).resid_qv)];
    else
      resid_qv_st=[];
    end
    
    if exist('supid','var')
      above_sup_st=num2str(round(nabove_sup/length(find(CL21.supdat(supid,:)==j))*1000)/10,'%3.1f/');
      above_sup_st=[':' num2str(round(nabove/size(CL21.dat,2)*1000)/10) '%-' above_sup_st(1:(end-1))];
    else
      above_sup_st=[];
    end
    
    if length(gidx)>1
      title(['P:' regexprep(st,'\(.*\)','') char(10) ... 
             'W:' regexprep(Wst,'\(.*\)','') char(10) ... 
             ' (' num2str(nabove) above_sup_st ')' char(10) ...
             'score:' num2str(CL21.ads(regs{k}(i).peak,k)) ' qv:' num2str(CL21.qv(regs{k}(i).peak,k)) ...
             resid_qv_st char(10) ...
             list2str({rgg(gidx).symbol},{','})],...
            'FontSize',6);
    elseif length(gidx)==1
      title(['P:' regexprep(st,'\(.*\)','') char(10) ... 
             'W:' regexprep(Wst,'\(.*\)','') char(10) ...
             ' (' num2str(nabove) above_sup_st ')' char(10) ...
             'score:' num2str(CL21.ads(regs{k}(i).peak,k)) ' qv:' num2str(CL21.qv(regs{k}(i).peak,k)) ...
             resid_qv_st char(10) ...
             rgg(gidx(1)).symbol ':' ...
             hum_num2str(double(rgg(gidx(1)).start),1) '-' hum_num2str(double(rgg(gidx(1)).end),1) ],...
            'FontSize',6);      
    else
      title(['P:' regexprep(st,'\(.*\)','') char(10) ... 
             'W:' regexprep(Wst,'\(.*\)','') char(10) ...
             ' (' num2str(nabove) above_sup_st ')' char(10) ...
             'score:' num2str(CL21.ads(regs{k}(i).peak,k)) ' qv:' num2str(CL21.qv(regs{k}(i).peak,k)) ...
             resid_qv_st],'FontSize',6);      
    end
    
    hold on
    if mod(i,5)==0
      print('-dpsc2',sfname,'-append');
      clf;
    end
  end
  if mod(i,5)~=0 
    print('-dpsc2',sfname,'-append');
    clf;
  end
  ps2pdf(sfname,['scatter_' kstr{k}  ext '.pdf'],0);
end
