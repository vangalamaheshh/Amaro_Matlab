function res=display_D_elem(gr,elem_type,Dord,gdend,sdend,varargin)
% display element in current axes

res=[];
switch elem_type
 case 'data'
  bluepink;
  if ~isempty(varargin)
    vtmp=varargin{1};
    Dtmp=preprocess_D(Dord,varargin{1});
    if ~isempty(Dtmp.dat)
      imagesc_trim(Dtmp.dat);
    end
  else
    if ~isempty(Dtmp.dat)
      imagesc_trim(Dord.dat);
    end
%    imagesc(Dord.dat);
%    imagesc_trim(dna_norm(Dord.dat')');    
  end
  ax=axis;
  % line([ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'Color','b','LineWidth',3);
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis on;
  box on;
  bluepink;
 
 case 'dataorig'
  % sample the copy number vertically
  maxypix = 10000; % maximum number of pixels
  samplint = ceil(size(Dord.dat,1)/maxypix);
  % draw heatmap
  imagesc(Dord.dat(1:samplint:end,:));
  if ~isempty(varargin)
    caxis(varargin{1});
  end
  % no axis ticks
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis on;
  box on;
  bluepink;
  if ~isempty(varargin)
    if length(varargin)>1
      colormap(varargin{2});
    end
  end
  res=caxis;
 
 case 'genotype'
  x1=ones(size(Dord.affy_calls));
  x2=ones(size(Dord.affy_calls));
  x3=ones(size(Dord.affy_calls));
  
  x1(Dord.affy_calls==1)=1;
  x2(Dord.affy_calls==1)=0;
  x3(Dord.affy_calls==1)=0;

  x1(Dord.affy_calls==2)=0;
  x2(Dord.affy_calls==2)=0;
  x3(Dord.affy_calls==2)=1;

  x1(Dord.affy_calls==3)=1;
  x2(Dord.affy_calls==3)=1;
  x3(Dord.affy_calls==3)=0;
  
  image(cat(3,x1,x2,x3));
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis on;
  box on;
  res=[];
  
 case 'consensusmat'
  imagesc(Dord.dat);
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis on;
  box on;
  whitered;

 case 'colormat'
  colormat(Dord.dat,varargin{:});
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis on;
  box on;

 case 'legendbar'
  draw_legendbar(varargin{:});
  
 case 'gacc'
  res=draw_names_box(strvcat(Dord.gacc),varargin{:});
  noticks;
  
 case 'gdesc'
  res=draw_names_box(strvcat(Dord.gdesc),varargin{:});
  noticks;
  
 case 'gdist'
  if isfield(gdend,'dist_mat')
    display_D_dist(gdend.dist_mat)
  else
    display_D_dist(Dord,'genes',gdend.params);
  end
  
 case 'gdend'
  if exist('gdend','var') && ~isempty(gdend)
    set(gca,'CameraUpVector',[-1 0 0]);
    set(gca,'YDir','reverse');
    if length(varargin)>0
      res=show_dend_lnk_idx(gdend.lnk,gdend.idx,[],ceil(0.01* ...
                                                        length(gdend.idx)),varargin{1},-ones(length(gdend.idx),1),varargin{2});
    else
      res=show_dend_lnk_idx(gdend.lnk,gdend.idx,[],ceil(0.01* ...
                                                        length(gdend.idx)),'reduced_big',-ones(length(gdend.idx),1));
    end    
    noticks;
  end
  
 case 'gsupdat'
  if isfield(Dord,'gsupdat')
    axis([0.5 size(Dord.gsupdat,2)+0.5 0 0.000001]);
    set(gca,'CameraUpVector',[-1 0 0]);
    for i=1:size(Dord.gsupdat,1)
      add_patches(Dord.gsupdat(i,:),Dord.gsupmark(i).colormap,0.8,Dord.gsupmark(i).height,0.5*(i>1));
    end
    noticks;
  end
  
 case 'gsupacc'
  % DID NOT TEST YET
  hvec=cat(1,Dord.gsupmark(:).height)+0.5;
  hvec(1)=hvec(1)-0.25;
  hvec(end)=hvec(end)-0.25;
  [suptitle,supnames]=break_sup_names(Dord.gsupacc);
  draw_names_box(strvcat(suptitle),varargin{1},varargin{2},hvec);
  noticks;
 
 
 case 'sdesc'
  res=draw_names_box(Dord.sdesc,varargin{:});
 
 case 'sdist'
  display_D_dist(Dord,'samples',sdend.params);
  
 case 'sdend'
  if exist('sdend','var') && ~isempty(sdend)
    if length(varargin)>0
      res=show_dend_lnk_idx(sdend.lnk,sdend.idx,[],1,varargin{1},- ...
                            ones(length(sdend.idx),1),varargin{2});
    else
      res=show_dend_lnk_idx(sdend.lnk,sdend.idx,[],1,'all',- ...
                            ones(length(sdend.idx),1));
    end    
    ax1=axis;
    axis tight
    ax2=axis;
    axis([0.5 length(sdend.idx)+0.5 ax2(3) ax1(4)]);
    noticks;
  end
  
 case 'ssupdat'
  if isfield(Dord,'supdat')
    if ~isfield(Dord,'supmark')
      c=read_colorscheme('~/CancerGenomeAnalysis/trunk/matlab/colorscheme.txt');
      Dord=add_supmark(Dord,c);
    end
    if nnz(Dord.supdat<0)
      disp('Fixing supdat');
      Dord=fix_supdat(Dord);
    end
    axis([0.5 size(Dord.supdat,2)+0.5 0 0.000001]);
    for i=1:size(Dord.supdat,1)
      if ~isfield(Dord.supmark(i),'patchwidth')
        pw=0.8;
      else
        pw=Dord.supmark(i).patchwidth;
      end
      if isfield(Dord.supmark(i),'linewidth')
        add_patches(Dord.supdat(i,:),Dord.supmark(i).colormap,pw, ...
                    Dord.supmark(i).height,0.5*(i>1),Dord.supmark(i).linewidth);
      else
        add_patches(Dord.supdat(i,:),Dord.supmark(i).colormap,pw, ...
                    Dord.supmark(i).height,0.5*(i>1));
      end
    end
    noticks;
  end
  
 case 'ssupacc'
  if isfield(Dord,'supdat')
    if ~isfield(Dord,'supmark')
      c=read_colorscheme('~/CancerGenomeAnalysis/trunk/matlab/colorscheme.txt');
      Dord=add_supmark(Dord,c);
    end
    hvec=cat(1,Dord.supmark(:).height)+0.5;
    hvec(1)=hvec(1)-0.25;
    hvec(end)=hvec(end)-0.25;
    [suptitle,supnames]=break_sup_names(Dord.supacc);
    if length(varargin)>2
      res=draw_names_box(strvcat(suptitle),varargin{1},varargin{2},hvec,varargin{3:end});
    else
      res=draw_names_box(strvcat(suptitle),varargin{1},varargin{2},hvec);
    end
  end
  noticks;
  
 case 'ssuplegend'
  if isfield(Dord,'supdat')
    if ~isfield(Dord,'supmark')
      c=read_colorscheme('~/CancerGenomeAnalysis/trunk/matlab/colorscheme.txt');
      Dord=add_supmark(Dord,c);
    end    
    draw_legend_box(Dord,varargin{:});
  end


 case 'colorbar'
  switch length(varargin)
   case 0
    draw_colorbar('hori');
   case 4
    draw_colorbar(varargin{1},varargin{2},varargin{3},varargin{4})
   case 6 % was 5 
    if strcmp(varargin{2},'for')
      r=get_subplotgrid_data_range(gr,varargin{3},varargin{4});
      draw_colorbar(varargin{1},[1+(63/varargin{5})*(0:varargin{5})],cellstr(num2str((r(1):(r(2)-r(1))/varargin{5}:r(2))'))',varargin{6})      
    else
      error('unknown colorbar param');
    end
   otherwise
      error('wrong number of parameters for colorbar');    
  end   
%  noticks;
  axis on
  box on 

 case 'scatterplot'
  scatter_plot(Dord,rc,supid1,supid2)
  axis on
  box on
  
 case 'text'
  res=text(varargin{:});
  noticks;  
  
 case 'chrn'
   if isfield(Dord,'chrn')
      [u,ui,uj]=unique_keepord(Dord.chrn);
      rl=runlength(uj);
      sz=rl(:,2)-rl(:,1)+1;
      res=draw_names_box(num2chromosome(unique_keepord(Dord.chrn)),varargin{1},varargin{2},sz,varargin{3:end});
   end 
   
 case 'cytoname'
    if isfield(Dord,'cyto_name') && isfield(Dord.cyton)
      %%% HERE add cyto band text
      [u,ui,uj]=unique_keepord(Dord.cyton);
      rl=runlength(uj);
      sz=rl(:,2)-rl(:,1)+1;
      res=draw_names_box(num2chromosome(unique_keepord(Dord.chrn)),varargin{1},varargin{2},sz,varargin{3:end});
   end 
 
 case 'chrcyto'
   if exist('varargin','var') && length(varargin)>0
     if isfield(Dord,'cyto_stain')
       % take care of NaNs
       tmp=Dord.cyto_stain;
       tmp(isnan(tmp))=0;
       tmp=tmp/100;
       image(repmat(1-tmp,[1 1 3]));
       noticks;
       axis on
       box on
     end
   elseif isfield(Dord,'chrn')
%       [Dord,supid]=add_D_sup(Dord,'chr','chr',mod(Dord.chrn',2),'rows');
%       Dord=reorder_D_sup(Dord,'rows',supid);
%       Dord.gsupmark(1).colormap=[1 1 1; 0 0 0];
%       Dord.gsupmark(1).height=1;
%       res=display_D_elem(gr,'gsupdat',Dord,gdend,sdend,varargin{:});
        image(repmat(mod(Dord.chrn,2)*0.9,[1 1 3]));
        if isfield(Dord,'armn')
          ax=axis;
          armpos=find(diff(Dord.armn)>0);
          for i=1:length(armpos)
            res(i)=line(ax(1:2),[armpos(i) armpos(i)]+0.5,'LineStyle','-','Color',[0.6 0.6 0.6],'LineWidth',1);
          end
        end
        noticks;            
    end
case 'fxascore'
    k=varargin{1};
    if iscell(Dord.fxa)
      Dord.fxa=cat(2,Dord.fxa{:});
    end
    if isfield(Dord,'fxa_extra')
      for ei=k:2:size(Dord.fxa_extra,2)
        ph=plot(Dord.fxa_extra(:,ei)); hold on
        set(ph,'Color',varargin{4}(floor((ei-k)/2)+1,:)); 
      end
    end
    plot(Dord.fxa(:,k)); hold on
    set(gca,'CameraUpVector',[-1 0 0]);

    ax=axis;
    if length(varargin)>1 && ~isempty(varargin{2})
      axis([0.5 length(Dord.fxa(:,k))+0.5 varargin{2}]);
    else
      axis([0.5 length(Dord.fxa(:,k))+0.5 ax(3:4)]);
    end
    if length(varargin)>2 && ~isempty(varargin{3})
      res(1)=line([1 length(Dord.fxa(:,k))],[ varargin{3} varargin{3}],'Color','g');
    end
    % add scale in each graph
    set(gca,'XTick',[]);
    yt=get(gca,'YTick');
    ytl=get(gca,'YTickLabel');
    text((0.5+length(Dord.fxa(:,k))*1.001)*ones(length(yt),1),yt,deblank(cellstr(ytl)),'HorizontalAlignment','center','VerticalAlignment','top');
    set(gca,'YTickLabel',[]);   
    box on;
    
 case 'reglines'
  ax=axis;
  rng=varargin{1};
  idx=varargin{2};
  reg=varargin{3};
  line(ax(1:2),reg.peak_st*ones(1,2)-rng(1)-0.5,'LineStyle','-','Color','k');
  line(ax(1:2),reg.peak_en*ones(1,2)-rng(1)+0.5,'LineStyle','-','Color','k');
  line(ax(1:2),reg.peak_wide_st*ones(1,2)-rng(1)-0.5,'LineStyle','--','Color','g');
  line(ax(1:2),reg.peak_wide_en*ones(1,2)-rng(1)+0.5,'LineStyle','--','Color','g');
  for i=1:length(idx)
    seg=reg.segments{idx(i)};
    if ~isempty(seg)
      line([ i-0.5 i+0.5 i+0.5 i-0.5 i-0.5],[seg.st-rng(1)+0.5 seg.st-rng(1)+0.5 seg.en-rng(1)+1.5 seg.en-rng(1)+1.5 ...
                          seg.st-rng(1)+0.5],'LineStyle','-','Color','k','LineWidth',0.1);
    end
  end
  
 case 'postick'
  plot(Dord.pos,Dord.pos);
  yt=get(gca,'YTick');
  ytl=get(gca,'YTickLabel');
  cla
  axis([0 1 0.5 size(Dord.dat,1)+0.5]);
  ax=axis;
  
 case 'pvs'
  if isfield(Dord,'pvs')
        k=varargin{1};
        use_loglog=varargin{2};
        q_thresh=varargin{3};
        if iscell(Dord.pvs)
          Dord.pvs=cat(2,Dord.pvs{:});
        end
        if use_loglog
          y=log10(-log10(Dord.pvs(:,k))+1);
        else
          y=-log10(q{k});
        end
        plot(y);
        set(gca,'CameraUpVector',[-1 0 0]);
      
        ax=axis;
        axis([1 length(Dord.pvs{k}) ax(3:4)]);
        if (0)
          yt=get(gca,'YTick');
          for i=1:length(yt)
            if ~isempty(Dord.fxa{k})
              tmp=Dord.fxa{k}(find(y<=yt(i)));
              if ~isempty(tmp)
                st(i)=max(tmp);
              else
                st(i)=NaN;
              end
              text(length(Dord.pvs{k}),yt(i),num2str(st(i)));
            end
          end
          
          if use_loglog
            set(gca,'YTick',yt,'YTickLabel',num2str(10.^(-((10.^(yt'))-1))));
          else
            set(gca,'YTick',yt,'YTickLabel',num2str(10.^(-(yt)')));
          end
        end
        
        if q_thresh>=0
          if use_loglog
            line([1 length(Dord.pvs{k})],[log10(-log10(q_thresh)+1) log10(-log10(q_thresh)+1)],'Color','r');
          else
            line([1 length(Dord.pvs{k})],[-log10(q_thresh) -log10(q_thresh)],'Color','r');
          end
        end
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        box on;
  end
  
end

%function noticks
%set(gca,'XTick',[]);
%set(gca,'YTick',[]);
%axis off
%box off







