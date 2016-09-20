function plot_region_regs(CL21,regs,ts,F1c,Fts,ssti,cyto,ext,rgg,use_regs_for_resid_qv,supid)
%% add support for supid and supmark
if ~exist('ext','var')
  ext=[];
end

if ~exist('use_regs_for_resid_qv','var')
  use_regs_for_resid_qv=0;
end


kstr={'amp','del'};

loc_chrn=cat(1,rgg.chrn);
loc_st=cat(1,rgg.start);
loc_en=cat(1,rgg.end);

score_type=struct('method','nxa','amp_thresh',ts(1),'del_thresh',-ts(2),'res',0.001);
[samp,sdel]=snp_score(CL21,score_type);
fxa_extra=[];
if exist('supid','var')
  [u,ui,uj]=unique(CL21.supdat(supid,:));
  if length(u)>1
    for uidx=1:length(u)
      [tmp1,tmp2]=snp_score(reorder_D_cols(CL21,find(uj==uidx)),score_type);
      fxa_extra=[fxa_extra tmp1 tmp2];
    end
  end
end
if ~isempty(fxa_extra)
  CL21.fxa_extra=fxa_extra;
end


disp_p=struct('x',struct('sizes',[0.3 0.1 0.03 0.03 0.1 3 0.5 ],'gaps',[1 1 0 1 1 0 1 1],'border',0.2),...
              'y',struct('sizes',[0.2 1.3 0.1],'gaps',[1 1 1 5],'border',0.2),...
              'items',...
              {{{1,1,'ssupacc','vert',8},...
                {2,1,'ssuplegend',2,8,6},...
                {1,6,'ssupdat'},...
                {2,6,'dataorig',[-1 1]},...
                {2,2,'chrn','vert',8,-0.5},...
                {2,3,'chrcyto'},...
                {2,4,'chrcyto',1},...
                {2,5,'postick'},...
                {3,6,'colorbar','horizontal','for',2,4,4,6},...
               }});

for k=1:2
  sfname=['region_' kstr{k}  ext '.ps' ];
  for i=1:length(regs{k})
    figure(1); clf;
%    show_samples=find((3-2*k)*CL21.dat(regs{k}(i).peak,:)>ts(k) | ...
%                      (3-2*k)*F1c.dat(regs{k}(i).peak,:)>Fts(k));
    show_samples=1:size(CL21.dat,2);
    X=reorder_D_cols(CL21,show_samples);
    X.fxa=[samp sdel];
    
    rng=[regs{k}(i).peak_wide_st regs{k}(i).peak_wide_en];
    rng_sz=rng(2)-rng(1)+1;
    rng=[floor(rng(1)-max(0.1*rng_sz,100)) ceil(rng(2)+max(0.1*rng_sz,100))];
    rng=[max(1,rng(1)) min(size(X.dat,1),rng(2))];
    Y=reorder_D_rows(X,rng(1):rng(2));
    
    display_D(Y,[],[],{disp_p,struct('items',{{{2,7,'fxascore',k,[],[],Y.supmark(supid).colormap},...
                        {2,6,'reglines',rng,ssti(show_samples),regs{k}(i),1}}})});
    print('-dpsc2',sfname,'-append');
    clf;
  end
  ps2pdf(sfname,['region_' kstr{k}  ext '.pdf'],0);
end

return

sfname=['region_' kstr{k} ext ];
  

for k=1:2
  sfname=['scatter_' kstr{k} ext '.ps'];
  if exist(sfname,'file')
    delete(sfname);  end
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
      title([regexprep(st,'\(.*\)','') char(10) ... 
             ' (' num2str(nabove) above_sup_st ')' char(10) ...
             'score:' num2str(CL21.ads(regs{k}(i).peak,k)) ' qv:' num2str(CL21.qv(regs{k}(i).peak,k)) ...
             resid_qv_st char(10) ...
             list2str({rgg(gidx).symbol},{','})],...
            'FontSize',6);
    elseif length(gidx)==1
      title([regexprep(st,'\(.*\)','') char(10) ...
             ' (' num2str(nabove) above_sup_st ')' char(10) ...
             'score:' num2str(CL21.ads(regs{k}(i).peak,k)) ' qv:' num2str(CL21.qv(regs{k}(i).peak,k)) ...
             resid_qv_st char(10) ...
             rgg(gidx(1)).symbol ':' ...
             hum_num2str(double(rgg(gidx(1)).start),1) '-' hum_num2str(double(rgg(gidx(1)).end),1) ],...
            'FontSize',6);      
    else
      title([regexprep(st,'\(.*\)','') char(10) ... 
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
