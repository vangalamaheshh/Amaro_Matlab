function [regs,pvs]=generate_regs(CL21,fsuffix,ts,h_amp,h_del,naamp,nadel,thst_amp,thst_del,pv_thresh,smooth_window_size,use_loglog);

if ~exist('smooth_window_size','var')
  smooth_window_size=501;
  disp(['Using smooth_window_size=' num2str(smooth_window_size)]);
end

if ~exist('use_loglog','var')
  use_loglog=1;
end

regs={};
pvs={};
% ts=[0.35 0.3];
% ts=[0.3];
les_acc={};
% ts=[1.16];
for k=1:length(ts)
  for dodel=0:1
    % dodel=1;

    if dodel
      thst=thst_del(:,k);
      name='deletions';
    else
      thst=thst_amp(:,k);
      name='amplifications';
    end
    if sum(thst)<=1.01 % it is a distribution
        chst=cumsum(thst)/sum(thst);
        one_chst=flipud(cumsum(flipud(thst)))/sum(thst);
    else
        chst=(cumsum(thst)+1)/(sum(thst)+2);
    end
    min(find(chst>0.95))
    semilogy(1-chst(1:200))
    
    if dodel
      y=nadel;
    else
      y=naamp;
    end
    
    y(isnan(y))=0;
    tail=1-[0; chst];
    tail=one_chst;
    s=y;
    s=floor(s*10);
    s(s>1000)=1000;
    pv=tail(s+1);

%    pv(pv<eps)=eps;
    
    semilogy(pv) 

    fdrpv=calc_fdr_value(pv);
    pvs{k,dodel+1}=fdrpv;
    
    reg=zeros(size(fdrpv,1),1);
    reg(fdrpv<pv_thresh)=1;
    
    % previously we were not dealing with the X well (prostates will be called deleted).
    % now we add 1 to all male samples so we do not need this
%    if dodel
%      reg(find(CL21.chrn==chromosome2num('X')))=0;
%    end
    
    %smoothreg=conv(reg,ones(200,1));
    if (0)
      CS=CL21; % Corig_filtered;
      CS.dat=reg;
      CS.smooth=[];
      CS.sdesc='x';
      smoothreg=smooth_copy_number(CS,smooth_window_size,0,'mean'); 
      smoothreg=smoothreg.smooth;
    else
      CS=CL21;
      CS.cbs=reg;
      CS.cbs_rl=[];
      CS.sdesc='x';
      CS.smooth=[];
      smoothreg=smooth_cbs(CS,smooth_window_size,1);
      smoothreg=smoothreg.cbs;
    end
    %smoothreg=reg;
    
    regions=[];
    for ci=1:max(CL21.chrn)
      inchr=find(CL21.chrn==ci);
      posreg=(smoothreg(inchr)>0);
      df=find(abs(diff([0; posreg; 0])));
      st=df(1:2:end);
      en=df(2:2:end)-1;
      regions=[regions; st+min(inchr)-1 en+min(inchr)-1 zeros(size(st,1),1)];
    end

    f=1;
    for i=1:size(regions,1)
      st=regions(i,1);
      en=regions(i,2);
      [mx,mi]=max(y(st:en));
      regions(i,3)=st+mi-1; % peak position
      pk=st+mi-1;
      pkrange=find(y(st:en)==mx);
      pkrange=pkrange+st-1;
      pkrl=runlength(double(y(st:en)==mx)');
      pkrl=pkrl(find(pkrl(:,3)==1),:);
      pkst=['['];
      for pi=1:size(pkrl,1)
        pkst=[pkst sprintf('%g-%g ',CL21.pos(pkrl(pi,1)+st-1),CL21.pos(pkrl(pi,2)+st-1))];
      end
      if size(pkrl,1)>1
        disp('WARNING: more than one segment');
      end
      regs{k,dodel+1}(i).peak=round(0.5*(pkrl(1,1)+st-1+pkrl(1,2)+st-1));
      regs{k,dodel+1}(i).peak_st=pkrl(1,1)+st-1;
      regs{k,dodel+1}(i).peak_en=pkrl(1,2)+st-1;      
      regs{k,dodel+1}(i).st=st;
      regs{k,dodel+1}(i).en=en;
      pkst=[pkst ']'];
      les_st=sprintf('%s:%g Mb - %s:%g Mb  peak at %s:%s, pv %f %f : ( %d:%d )',...
                 CL21.chr{st},CL21.pos(st),CL21.chr{en},CL21.pos(en),CL21.chr{pk}, ...
                 pkst,pv(pk),fdrpv(pk),st,en); %CL21.pos(pk)
      fprintf(f,les_st);
      les_acc{i}=les_st;
      fprintf(f,'\r\n');
    end
    
    for i=1:max(CL21.chrn)
      chrnpos(i)=round(mean(find(CL21.chrn==i)));
    end

    figure(1); clf;
    subplot(6,1,1:3);
    if use_loglog
      plot(log10(-log10(fdrpv)+1));
    else
      plot(-log10(fdrpv));
    end
    ax=axis;
    axis([1 length(fdrpv) ax(3:4)]);
    yt=get(gca,'YTick');
    if use_loglog
      set(gca,'YTick',yt,'YTickLabel',num2str(10.^(-((10.^(yt'))- ...
                                                     1))));
    else
      set(gca,'YTick',yt,'YTickLabel',num2str(10.^(-(yt)')));
    end
    set(gca,'XTick',chrnpos,'XTickLabel',CL21.chr(chrnpos));
    axs(1)=gca;
    if use_loglog
      line([1 length(fdrpv)],[log10(-log10(pv_thresh)+1) log10(-log10(pv_thresh)+1)],'Color','r');
    else
      line([1 length(fdrpv)],[-log10(pv_thresh) -log10(pv_thresh)],'Color','r');
    end
    subplot(6,1,4);
    imagesc(reg'>0);
    axs(2)=gca;
    set(gca,'YTick',[]);
    set(gca,'XTick',chrnpos,'XTickLabel',CL21.chr(chrnpos));
    subplot(6,1,5);
    colormap gray
    imagesc(smoothreg'>0);
    colormap gray
    axs(3)=gca;
    set(gca,'YTick',[]);
    set(gca,'XTick',chrnpos,'XTickLabel',CL21.chr(chrnpos));
    subplot(6,1,6);
    colormap gray
    imagesc(mod(CL21.chrn,2)');
    axs(4)=gca;
    set(gca,'YTick',[]);
    set(gca,'XTick',chrnpos,'XTickLabel',CL21.chr(chrnpos));
    linkaxes(axs,'x');
    
    print_pdf([ name '_'  num2str(size(CL21.dat,2)) '_' num2str(ts(k)) ...
                '_' fsuffix '_figure']);
    %t2=log2(3)-1;
    if dodel
      %  t1=-0.6; % log2(4/3)-1
      t1=-ts(k); % log2(4/3)-1
      write_eisen_dat(['deletions_log_' num2str(size(CL21.dat,2)) '_' num2str(ts(k)) '_' fsuffix '.txt'],strvcat(les_acc),repmat(' ',length(regions),1),strvcat(CL21.sdesc),'Deletions', ...
                      double(CL21.smooth(regions(:,3),:)<t1),[],[],1);
    else
      %  t1=0.6; % log2(3)-1
      t1=ts(k); % log2(3)-1
      write_eisen_dat(['amplifications_log_' num2str(size(CL21.dat,2)) '_' num2str(ts(k)) '_' fsuffix '.txt'],strvcat(les_acc),repmat(' ',length(regions),1),strvcat(CL21.sdesc),'Amplifications', ...
                      double(CL21.smooth(regions(:,3),:)>t1),[],[],1);
    end
  end
end
