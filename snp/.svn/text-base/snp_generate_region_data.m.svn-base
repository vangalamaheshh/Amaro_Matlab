function tab=snp_generate_region_data(CL21r,A,u95ll,rg,cyto,regs,calls,ts,pvs,MD,ext)

if ~exist('ext','var') 
  ext=[];
end

if mean(CL21r.pos)<10000 % is in Mb
  error('DOES NOT WORK WITH MB ANYMORE');
end

if ~isempty(A) && ischar(A.gacc)
  A.gacc=cellstr(A.gacc);
end

tab={};
for k=1:length(calls)
  for pki=1:size(calls{k},1)  % :-1:1
    disp([ k pki ]);
    tab{k}.idx(pki)=pki;
    curc=calls{k}(pki,:);
    cls1=find(curc==1);
    cls0=find(curc==0);
    
    if length(cls0)==0 || length(cls1)==0 
      disp('One of the classes is empty');
      continue;
    end
    
    if isfield(regs{k}(pki),'peak_wide_st')
      regs{k}(pki).enl_peak_st=regs{k}(pki).peak_wide_st;
      regs{k}(pki).enl_peak_en=regs{k}(pki).peak_wide_en;
    else
      [regs{k}(pki).enl_peak_st,regs{k}(pki).enl_peak_en]=...
          enlarge_range(CL21r,regs{k}(pki).peak_st,regs{k}(pki).peak_en,15,15,...
                        regs{k}(pki).st,regs{k}(pki).en);
    end
    
    [tmp1,tmp2,enl_region_to_next_snp]=genomic_location(CL21r,{[regs{k}(pki).enl_peak_st ...
                        regs{k}(pki).enl_peak_en]},cyto,1);
    [tmp1,tmp2,peak_region_to_next_snp]=genomic_location(CL21r,{[regs{k}(pki).peak_st ...
                        regs{k}(pki).peak_en]},cyto,1);
    
    genes_in_region=genes_at(rg,CL21r.chrn(regs{k}(pki).enl_peak_st), ...
                             enl_region_to_next_snp{1}(1), enl_region_to_next_snp{1}(2));
    genes_in_peak=genes_at(rg,CL21r.chrn(regs{k}(pki).peak_st), ...
                           peak_region_to_next_snp{1}(1),peak_region_to_next_snp{1}(2));
    tab{k}.chr(pki)=CL21r.chrn(regs{k}(pki).peak_st);
    tab{k}.region_start(pki)=CL21r.pos(regs{k}(pki).st);
    tab{k}.region_end(pki)=CL21r.pos(regs{k}(pki).en);
    
    tab{k}.peak_start(pki)=CL21r.pos(regs{k}(pki).peak_st);
    tab{k}.peak_end(pki)=CL21r.pos(regs{k}(pki).peak_en);
    
    tab{k}.enlarged_peak_start(pki)=CL21r.pos(regs{k}(pki).enl_peak_st);
    tab{k}.enlarged_peak_end(pki)=CL21r.pos(regs{k}(pki).enl_peak_en);
    
    tab{k}.n_genes_in_region(pki)=length(genes_in_region);
    tab{k}.n_genes_in_peak(pki)=length(genes_in_peak);

    if ~isempty(genes_in_region)
      tab{k}.genes_in_region{pki}=sprintf('%s,',rg(genes_in_region).symbol);
    else
      tab{k}.genes_in_region{pki}='Empty';
    end
    
    if ~isempty(genes_in_peak)
      tab{k}.genes_in_peak{pki}=sprintf('%s,',rg(genes_in_peak).symbol);
    else
      tab{k}.genes_in_peak{pki}='Empty';
    end
    
    currg={rg(genes_in_region).locus_id};  
    if ~isempty(currg) && ~ischar(currg{1})
      currg=cellstr(num2str(cat(1,currg{:}),'%d'));
    end
    
    curA=[];
    curAord=[];
    P=[];

    if ~isempty(A)
      if exist('u95ll_h','var')
        disp('Here');
        [M,m2,m1]=match_string_sets_hash(u95ll,currg,u95ll_h,u95ll_us1j);
      else
        [M,m2,m1,u95ll_h,u95ll_us1j]=match_string_sets_hash(u95ll,currg);
      end
      [genes_on_chip,ui,uj]=unique(m2);
      m1=m1(ui);
      cur_genes=rg(genes_in_region(m1));
      
      tab{k}.n_genes_on_chip(pki)=length(cur_genes);
      if ~isempty(cur_genes)
        tab{k}.genes_on_chip{pki}=sprintf('%s,',cur_genes(:).symbol);
      else
        tab{k}.genes_on_chip{pki}='Empty';
      end
      
      not_on_chip=rg(genes_in_region(setdiff(1:length(genes_in_region),m1)));

      curA=reorder_D_rows(A,genes_on_chip);
      curA.gdesc=remove_quotes({cur_genes(:).gene});
      curA.gsymb={cur_genes(:).symbol};
      curA.grg=cur_genes;
      curA.chrn=cat(1,cur_genes(:).chrn);
      curA.start=cat(1,cur_genes(:).start);
      curA.end=cat(1,cur_genes(:).end);
      mid=0.5*(curA.start+curA.end);
      [mids,midi]=sort(mid);
      curA=reorder_D_rows(curA,midi);
      
      % should be 2-sided and then pick the correct direction
      % correct for multiple hypotheses
      if (k==1) && (length(calls)==2)
        [P,S]=differential_analysis(curA,cls1,cls0,struct('method','ttest1side_minvar','minvar',0.5^2),0);
      else
        [P,S]=differential_analysis(curA,cls0,cls1,struct('method','ttest1side_minvar','minvar',0.5^2),0);
      end
      % to take care of cases with 1 samples vs the rest
      
      P(isnan(P))=0.1;
      [ss,si]=sort(S);
      [sp,spi]=sort(P);
      if ~isempty(P)
        xx=curA.gsymb(spi(1:min(length(spi),3)));
        tab{k}.top3{pki}=sprintf('%s,',xx{:});
      else
        tab{k}.top3{pki}='Empty';
      end

      % More Data 
      if exist('MD','var') && ~isempty(MD)
        for mdi=1:length(MD)
          [Px,Sx]=differential_analysis(MD{mdi},cls1,cls0,struct('method','ttest'),0);
          Px(isnan(Px))=0.1;
          [ssx,six]=sort(Sx);
          [spx,spix]=sort(Px);
          if ~isempty(Px)
            xx={};
            for jj=1:min(15,length(spix))
              xx{jj}=[ deblank(MD{mdi}.gdesc(spix(jj),:)) ' [' ...
                       num2str(Sx(spix(jj))) ',' num2str(Px(spix(jj))) ...
                       ']' ];
            end
            tab{k}.md{pki,mdi}=sprintf('%s,',xx{:});
            tab{k}.mdres{pki,mdi}={Sx,Px};
          else
            tab{k}.md{pki,mdi}='Empty';
            tab{k}.mdres{pki,mdi}={};           
          end              
        end
      end
      
    else
      cls0=1;
      cls1=2:size(CL21r.dat,2);
      not_on_chip=rg(genes_in_region);
      tab{k}.n_genes_on_chip(pki)=0;
      tab{k}.genes_on_chip{pki}='Empty';
      tab{k}.top3{pki}='Empty';
      tab{k}.md={};
      genes_on_chip=[];
    end
    
    if 1 || (~isempty(curA.dat) & 1)
      figure(1);clf;
      set(gcf,'Visible','off');
      if ~isempty(curA)
        curAord=reorder_D_cols(curA,[cls0 cls1]);
        X=curAord;
        v=zeros(1,size(X.dat,2));
        v((end-length(cls1)+1):end)=1;
        X=add_D_sup(X,num2str(pki),num2str(pki),v,'cols');
        c=read_colorscheme('~/projects/miRNA/data/colorscheme.txt');
        X=add_supmark(X,c);
        X.sdesc=X.sdesc(:,1:min(10,size(X.sdesc,2)));
        X.gacc=X.gsymb;
        [res,gr]=display_D(X,[],[],'datsupnames',[0 0.4 1 0.6]);
      end
      if exist('gr','var')
        pos=gr{4,3}.position;
        subplot('position',[pos(1) 0.01 pos(3) pos(2)-0.02]);
      else
        subplot(5,1,4:5);
      end
      
      imagesc_trim(CL21r.dat(regs{k}(pki).st:regs{k}(pki).en,[cls0 cls1]));
      bluepink;

      ns=size(CL21r.dat,2);
      ys=[regs{k}(pki).peak_st-regs{k}(pki).st+0.5 ...
          regs{k}(pki).peak_en-regs{k}(pki).st+1.5 ];
      ys2=[regs{k}(pki).enl_peak_st-regs{k}(pki).st+0.5 ...
           regs{k}(pki).enl_peak_en-regs{k}(pki).st+1.5 ];
      t=0.02;
      lw=1;
      line([t ns-t ns-t t t]+0.5,[ys(1) ys(1) ys(2) ys(2) ys(1)],'Color','black','LineWidth',lw);
      line([t ns-t ns-t t t]+0.5,[ys2(1) ys2(1) ys2(2) ys2(2) ys2(1)],'Color','Yellow','LineWidth',lw);
      disp([pki length(currg) length(genes_on_chip)]);
      print_D(['peak_' num2str(k) '_' num2str(pki) ext],{{'png','-r180'}},0);

      
      figure(2); clf;
      set(gcf,'Visible','off');
      x_peak_range=[ tab{k}.peak_start(pki) tab{k}.peak_end(pki)];
      x_range=[ tab{k}.enlarged_peak_start(pki) tab{k}.enlarged_peak_end(pki)];
      
      if isempty(P)
        top=0.3;
      else
        top=-log10(min(P))+0.3;
      end
      axis([ x_range+[-1 1] -0.3 top]);
      set(gca,'CameraUpVector',[-1 0 0]);
      
      ax=axis;
      %    patch([ x_peak_range x_peak_range(2:-1:1) x_peak_range(1) ],[ax(3) ...
      %                        ax(3) ax(4) ax(4) ax(3)],[0.8 0.8 ...
      %                        0.8],'FaceAlpha',0.8);

%    tab{k}.chr(pki)=CL21r.chrn(regs{k}(pki).peak_st);
%    tab{k}.region_start(pki)=CL21r.pos(regs{k}(pki).st);
%    tab{k}.region_end(pki)=CL21r.pos(regs{k}(pki).en);
    
%    tab{k}.peak_start(pki)=CL21r.pos(regs{k}(pki).peak_st);
%    tab{k}.peak_end(pki)=CL21r.pos(regs{k}(pki).peak_en);
    
%    tab{k}.enlarged_peak_start(pki)=CL21r.pos(regs{k}(pki).enl_peak_st);
%    tab{k}.enlarged_peak_end(pki)=CL21r.pos(regs{k}(pki).enl_peak_en);
       
      yc=CL21r.pos(regs{k}(pki).enl_peak_st:regs{k}(pki).enl_peak_en);
      %max((regs{k}(pki).peak_st-15),1):(min(regs{k}(pki).peak_en+ ...
      %                                               15,length(CL21r.pos))));
      xc=[-0.3 top];
      xx=CL21r.dat(regs{k}(pki).enl_peak_st:regs{k}(pki).enl_peak_en,:);
      %max((regs{k}(pki).peak_st-15),1):min((regs{k}(pki).peak_en+ ...
      %                                               15),length(CL21r.pos)),:);
      if (k==1) 
        xx(xx<ts(1))=0;
      else
        xx(xx>-ts(2))=0;
        xx=-xx;
      end
      sxx=sum(xx,2);
      sxx=-log10(pvs{k}(regs{k}(pki).enl_peak_st:regs{k}(pki).enl_peak_en));
      %max((regs{k}(pki).peak_st-15),1):min((regs{k}(pki).peak_en+15),length(CL21r.pos))));
      if length(yc)>1
        PH=pcolor(yc,xc,repmat(sxx,1,2)');
        set(PH,'FaceAlpha',0.8);
        cm=1-repmat([0:0.4/63:0.4],3,1)';
        if (k==1)  % && (length(calls)==2)
          cm(:,1)=1;
        else
          cm(:,3)=1;
        end         
        colormap(cm);
        caxis([0 max(sxx)]);
        shading('interp');
      end
      line([ x_peak_range([1 1]) ],[ax(3:4)],'Color','black');
      line([ x_peak_range([2 2]) ],[ax(3:4)],'Color','black');
      
      
      if ~isempty(curAord)
        for j=1:size(curAord.dat,1)
          line([curAord.grg(j).start curAord.grg(j).end],[-log10(P(j)) -log10(P(j)) ]);
          text(double(clip_to_range(0.5*(curAord.grg(j).start+curAord.grg(j).end),x_range)),...
               -log10(P(j))+top*0.01,deblank(curAord.gsymb{j}),'HorizontalAlignment','left');
          disp([ num2str(j) ') ' curAord.gacc{j} ' ' curAord.gsymb{j} ...
                 ' ' curAord.gdesc{j} ' ' num2str(-log10(P(j)))]);
        end
      end
      for j=1:length(not_on_chip)
        y1=-j/length(not_on_chip)*0.3;
        line([not_on_chip(j).start not_on_chip(j).end],[y1 y1 ],'Color','red');
        text(double(clip_to_range(0.5*(not_on_chip(j).start+not_on_chip(j).end),x_range)), ...
             y1+top*0.01,deblank(not_on_chip(j).symbol),'HorizontalAlignment','left');
      end
      set(gca,'CameraUpVector',[-1 0 0]);
      ypos=get(gca,'YTick');
      set(gca,'YTickLabel',[]);
      ax=axis;
      for i=1:length(ypos)
        th(i)=text(ax(2)+(ax(2)-ax(1))*0.02,ypos(i),num2str(ypos(i)),'HorizontalAlignment','center');
      end
      
      print_D(['genes_' num2str(k) '_' num2str(pki) ext],{{'png','-r180'}},0);

      figure(3); clf;
      set(gcf,'Visible','off');
      rng=regs{k}(pki).enl_peak_st:regs{k}(pki).enl_peak_en;
      %max((regs{k}(pki).peak_st-15),1):min((regs{k}(pki).peak_en+15),length(CL21r.pos));
      subplot(3,1,1);
      % raw
      %    xx=CL21r.raw(rng,:);
      %    imagesc_trim(xx);
      %    ca=caxis;
      %    colorbar
      %    bluepink;
      if isfield(CL21r,'raw')
        subplot(3,1,2);
        % batch corrected
        xx=CL21r.raw(rng,:);
        if abs(mean(CL21r.raw(:)))>0.5
          xx(xx<0.1)=0.1;
          xx=log2(xx)-1;
        end
        imagesc_trim(xx);
        ca=caxis;
        colorbar
        bluepink;
      end
      % smoothed
      if isfield(CL21r,'sm3')
        subplot(3,1,3);
        xx=CL21r.sm3(rng,:);
        imagesc_trim(xx);
        caxis(ca);
        colorbar
        bluepink;
        print_D(['smoothed_' num2str(k) '_' num2str(pki) ext],{{'png','-r180'}},0);
      end
      
      if exist('MD','var') && ~isempty(MD)
        for mdi=1:length(MD)
          figure(3+mdi); clf;
          set(gcf,'Visible','off');
           Sx=tab{k}.mdres{pki,mdi}{1};
           Px=tab{k}.mdres{pki,mdi}{2};
           [spx,spix]=sort(Px);
           l=min(15,length(spix));
           [ssx,six]=sort(Sx(spix(1:l)));
           iii=spix(six);
           X=reorder_D_rows(MD{mdi},iii);
           X=reorder_D_cols(X,[cls0 cls1]);
           v=zeros(1,size(X.dat,2));
           v((end-length(cls1)+1):end)=1;
           X=add_D_sup(X,num2str(pki),num2str(pki),v,'cols');
           c=read_colorscheme('~/projects/miRNA/data/colorscheme.txt');
           X=add_supmark(X,c);
           % X.sdesc=X.sdesc(:,1:10);
           display_D(X,[],[],'datsupnames');

           print_D(['markers_' num2str(k) '_' num2str(pki) '_' num2str(mdi) ext],{{'png','-r180'}},0);
        end
      end
    end 
    
    %  dist_from_peak=CL21r.pos(regs{k}(pki).peak)*1e6+5-...
    %      (cat(1,curAord.grg(:).start)+cat(1,curAord.grg(:).end))/2;
  end
end

for k=1:length(calls)   
  f=fopen(['table_' num2str(k) ext '.txt'],'w');
  fprintf(f,['index\tchromosome\tregion_start\tregion_end\tpeak_start\tpeak_end\tenlarged_peak_start\tenlarged_peak_end\t' ...
             'n_genes_in_region\tgenes_in_region\tn_genes_in_peak\tgenes_in_peak\tn_genes_on_chip\' ...
             'tgenes_on_chip\ttop 3\n']);
  if isfield(tab{k},'idx')
    for i=1:length(tab{k}.idx)
      if ~isfield(tab{k},'chr')
        fprintf(f,'%d',tab{k}.idx);
      else
        fprintf(f,['%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t' ...
                   '%d\t%s\t%d\t%s\t%d\t%s\t%s'],...
                tab{k}.idx(i),...
                num2chromosome(tab{k}.chr(i)),...
                tab{k}.region_start(i),...
                tab{k}.region_end(i),...
                tab{k}.peak_start(i),...
                tab{k}.peak_end(i),...
                tab{k}.enlarged_peak_start(i),...
                tab{k}.enlarged_peak_end(i),...
                tab{k}.n_genes_in_region(i),...
                tab{k}.genes_in_region{i},...
                tab{k}.n_genes_in_peak(i),...
                tab{k}.genes_in_peak{i},...                
                tab{k}.n_genes_on_chip(i),...
                tab{k}.genes_on_chip{i},...
                tab{k}.top3{i});
      end
      if isfield(tab{k},'md')
        for jj=1:size(tab{k}.md,2)
          fprintf(f,'\t%s',tab{k}.md{i,jj});
        end
      end
      fprintf(f,'\n');
    end
  end
  fclose(f);
end


