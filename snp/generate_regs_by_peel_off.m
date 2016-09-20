function [regs,C2]=generate_regs_by_peel_off(C,ads,d,q,score_type,score_thresh,wide_type,iter_type,peeloff_type)
%function call: [regs,C2]=generate_regs_by_peel_off(C,ads,d,q,score_type,score_thresh,wide_type,iter_type)
%old function call was: [regs,C2] =
%generate_regs_by_peel_off(C,d,q,score_type,score_thresh_method,score_type_smooth_
%sz,wide_type)  (I think!) --- code has been hopefully modified for
%backward compatibility, although change old code s.t. arg 2 has ads --
%see after_segmentation)
%---
% $Id$
% $Date: 2007-09-18 13:20:16 -0400 (Tue, 18 Sep 2007) $
% $LastChangedBy: rameen $
% $Rev$

if exist('wide_type','var') && isnumeric(wide_type)
    smooth_sz = wide_type;
end
if ~exist('wide_type','var') || isempty(wide_type)
  wide_type=struct('method','origloo');
elseif ischar(wide_type)
  tmp.method=wide_type;
  wide_type=tmp;
elseif isnumeric(wide_type)  % added for backward compatibility
    wide_type = iter_type;
end

if ~exist('iter_type','var') || isempty(iter_type) || exist('smooth_sz','var') %last bool added for backward compatibility
  iter_type=struct('method','no_iter');
elseif ischar(iter_type)
  tmp.method=iter_type;
  iter_type=tmp;
end

if ~exist('peeloff_type','var') || isempty(peeloff_type)
  peeloff_type=struct('method','segment');
end


if isnumeric(score_type)  %fix score_type for backwards compatibility
    tmp.method = 'nxa';
    tmp.amp_thresh = score_type(1);
    tmp.del_thresh = score_type(2);
    score_type = tmp;
    if isnumeric(wide_type)  % wide_type argument used to be smooth_sz in prior version (????)
        score_type.smooth_sz = smooth_sz;
    end
end

  



verbose('Using Peel Off',10);

show_plots=0;

% [ads{1},ads{2}]=snp_score(C,score_type);

regs=cell(1,2);
if ischar(C.sdesc)
  C.sdesc=cellstr(C.sdesc);
end

C2=C;

if strcmp('robust',wide_type.method)
  ts = [score_type.amp_thresh -1*score_type.del_thresh];
  conf_level = wide_type.conf_level;
  if isfield(wide_type,'ampdel')
    ampdel = wide_type.ampdel;
  else
    ampdel = [1 1];
  end
  perm_ads = generate_permuted_scores(C.dat,ts,wide_type.nperm,ampdel);
  range_dist = cell(1,2);
end

switch iter_type.method
  case 'no_iter'
   for k=1:2
     if k==1
       verbose('Amp',10);
     else
       verbose('Del',10);
     end
     range_dist{k} = cell(1,length(find(C.chrn == 1)));
     for ch=1:max(C.chrn)
       verbose(['Chr ' num2chromosome(ch)],20);
       
       in_chr=find(C.chrn==ch);
       if isempty(in_chr)
         continue;
       end
       chr_zero=min(in_chr)-1;
       chr_max = max(in_chr); % Marks last snp on current chromosome
       
       switch score_type.method
        case 'nxa'
         x=C.dat(in_chr,:);
         if k==1 % amp
           th=score_type.amp_thresh;
         else % del
           x=-x;
           th=-score_type.del_thresh;
         end
         x(x<th)=0;
        otherwise
         error('No such score method');
       end
       
       if nnz(isnan(x))>0
         tmp=make_D(x);
         tmp.chrn=C.chrn(in_chr);
         tmp.pos=C.pos(in_chr);
         tmp=impute_missing_values(tmp,'max');
         if isfield(score_type,'smooth_sz')
           tmp.cbs=tmp.dat;
           tmp=smooth_cbs(tmp,score_type.smooth_sz);
           tmp.dat=tmp.cbs;
         end
         x=tmp.dat;
       end

       y=(x>th);
       sc=nanmean(x,2);
       orig_sc=sc;

       rl=runlength(y);
       ns=size(y,2);

       % define mx_delta
       usc=unique(x(:));
       if length(usc)>1
         mx_delta=min(0.01/ns,min(0.5*diff(usc)));
       else
         mx_delta=0.01/ns;
       end
       verbose(['max delta = ' num2str(mx_delta)],30);
       
       [mx,mi]=max(sc);       
       verbose(['max(score)=' num2str(mx) '(?> ' num2str(score_thresh(k)) ')' ],20);
       
       if (show_plots==1)
         figure(1); clf;
         subplot(2,1,1);
         plot(sc); hold on;
         ax=axis;
       end
       
       [sads,sadsi]=sort(ads{k});
       sq=log(q{k}(sadsi));
       
       if show_plots==2
         X=reorder_D_rows(C2,in_chr);
         X.dat=x*[1-2*(k==2)];
         X.fxa{1}=sc;
         X.fxa{2}=sc;
         display_D(X,[],[],{'snp_fxa',struct('items',{{{2,5,'fxascore',1,[0 max(orig_sc)],score_thresh(k)}}})});
       end
       
       while mx>=score_thresh(k)
         n=length(regs{k})+1;
         
         % find samples 
         % HANDLE PEAKS WITH MORE THAN ONE SEGMENT!!!
         rg=find(sc>=mx-mx_delta);
         if any(diff(rg)>1)
           warning('the peak has more than one segment. should fix enlarged peak');
         end
         samples=find(y(mi,:));
         
         %      CC.sdesc=C.sdesc; CC.dat=C2.dat;
         %      save CC.x.mat CC
         
         % find segments in samples
         switch peeloff_type.method
          case 'segment'
           rli=find_rl(rl(samples),mi);
          case 'regions'
           rli=find_rl(rl(samples),mi); % find only overlapping segment
           regions_in_ch=find(peeloff_type.regions_rl(:,4)==ch);
           regions_rli=regions_in_ch(find_rl(peeloff_type.regions_rl(regions_in_ch,:),mi+chr_zero));
         end
           
         regs{k}(n).chrn=ch;
         v=zeros(1,ns);
         v(samples)=1;
         regs{k}(n).samples=v;
         
         
         regs{k}(n).peak_st=min(rg)+chr_zero;
         regs{k}(n).peak_en=max(rg)+chr_zero;
         regs{k}(n).peak=round(0.5*(regs{k}(n).peak_st+regs{k}(n).peak_en));
         regs{k}(n).qv=q{k}(in_chr(mi));
         regs{k}(n).resid_qv=exp(interp_pwl(sads,sq,mx));
         regs{k}(n).st=max([0; find(orig_sc(1:mi)<score_thresh(k)) ])+1+chr_zero;
         regs{k}(n).en=min([ find(orig_sc(mi:end)<score_thresh(k)); length(sc)-mi+1 ])-1+mi+chr_zero-1; 
         %% added -1 at end of prev. line since we want the last one that is still above score
         regs{k}(n).score=mx;
         verbose(['peak at ' genomic_location(C,{regs{k}(n).peak_st:regs{k}(n).peak_en})],20);
         
         z=zeros(size(x));
         for i=1:length(samples)
           for j=1:length(rli{i})
             rng=rl{samples(i)}(rli{i}(j),1):rl{samples(i)}(rli{i}(j),2);
             z(rng,samples(i))=x(rng,samples(i));
           end
         end
         %      if (k==2) && (ch==14)
         %        keyboard
         %      end
         
         % commented since we want to use the original wide region
         % use next line if you want to use all the data to look for the next boundary
         switch wide_type.method
          case 'robust'
           regs{k}(n).robust_reg_st = max([1; find(orig_sc(1:mi)<(mx- ...
                                                             score_thresh(k)))])+chr_zero;
           right_half = (mi:length(sc))+1;
           regs{k}(n).robust_reg_en = min([ right_half(find(orig_sc(mi:end)< ...
                                                            (mx-score_thresh(k))))+chr_zero-1 max(in_chr)]); 
           [wide_start wide_end thresh range_dist{k}] = robust_peak(regs{k}(n),ads{k},perm_ads{k},mx-score_thresh(k),score_type,regs{k}(n).robust_reg_st,regs{k}(n).robust_reg_en,conf_level(k),range_dist{k});
           
           regs{k}(n).peak_wide_st = max(chr_zero+1,wide_start-1);
           regs{k}(n).peak_wide_en = min(chr_max,wide_end+1);
          case 'origloo'
           [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en]=add_wide_peak(z,chr_zero);
           %         [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en]=add_wide_peak(regs{k}(n),x,chr_zero);
          case 'leave-k-out'
           %         if (k==2)
           %           keyboard
           %         end
           [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en]=add_wide_peak_k(z,chr_zero,wide_type.k);
          otherwise 
           error('no such method');
         end
         [regs{k}(n).peak_wide_st regs{k}(n).peak_st regs{k}(n).peak_en regs{k}(n).peak_wide_en];
         
         switch peeloff_type.method
           case 'segment'
            regs{k}(n).segments=cell(1,ns);
            for i=1:length(samples)
              regs{k}(n).segments{samples(i)}.rli=rli{i};
              regs{k}(n).segments{samples(i)}.st=rl{samples(i)}(rli{i},1)+chr_zero;
              regs{k}(n).segments{samples(i)}.en=rl{samples(i)}(rli{i},2)+chr_zero;
              pos=rl{samples(i)}(rli{i},1):rl{samples(i)}(rli{i},2);
              regs{k}(n).segments{samples(i)}.vec=x(pos,samples(i));
              
              verbose(['removing segment from ' C.sdesc{samples(i)} ' at ' ...
                       genomic_location(C,{regs{k}(n).segments{samples(i)}.st:regs{k}(n).segments{samples(i)}.en})],20);
              
              % replace with normal
              rl{samples(i)}=replace_rl(rl{samples(i)},rli{i},0); 
              x(pos,samples(i))=0;        
              y(pos,samples(i))=0;
              C2.dat(in_chr(pos),samples(i))=0;
              
              sc(pos)=sc(pos)-regs{k}(n).segments{samples(i)}.vec./ns;
              
              if show_plots>0
                subplot(2,1,1);
                plot(sc); hold on;
                subplot(2,1,2);
                tmp=zeros(length(sc),1);
                tmp(pos)=regs{k}(n).segments{samples(i)}.vec./ns;
                if show_plots==1
                  plot(tmp,'r-');
                  axis(ax);
                  disp([ 'min=' num2str(min(sc)) '; mean=' num2str(mean(regs{k}(n).segments{samples(i)}.vec./ns))]);
                  pause;
                end
              end
              if show_plots==2
                X=reorder_D_rows(C2,in_chr);
                X.dat=x*[1-2*(k==2)];
                X.fxa{1}=sc;
                X.fxa{2}=sc;
                display_D(X,[],[],{'snp_fxa',struct('items',{{{2,5,'fxascore',1,[0 max(orig_sc)],score_thresh(k)}}})});
              end
            end
           case 'regions'
            regs{k}(n).segments=cell(1,ns);
            for i=1:length(samples)
              regs{k}(n).segments{samples(i)}.rli=rli{i};
              rli_removed=find_rl(rl{samples(i)},peeloff_type.regions_rl(regions_rli,1:2)-chr_zero,0);
              if length(rli_removed)==0
                verbose(['Already removed segment in ' C.sdesc{samples(i)}],20);
                regs{k}(n).segments{samples(i)}.rli_removed=[];
                regs{k}(n).segments{samples(i)}.st=NaN;
                regs{k}(n).segments{samples(i)}.en=NaN;
                regs{k}(n).segments{samples(i)}.vec=[];
                continue;
              end
              regs{k}(n).segments{samples(i)}.rli_removed=rli_removed;
              st=min(rl{samples(i)}(rli_removed,1));
              en=max(rl{samples(i)}(rli_removed,2));
              regs{k}(n).segments{samples(i)}.st=st+chr_zero;
              regs{k}(n).segments{samples(i)}.en=en+chr_zero;
              pos=st:en;
              regs{k}(n).segments{samples(i)}.vec=x(pos,samples(i));
              
              verbose(['removing ' num2str(length(rli_removed)) ' segment(s) from ' C.sdesc{samples(i)} ' at ' ...
                       genomic_location(C,{regs{k}(n).segments{samples(i)}.st:regs{k}(n).segments{samples(i)}.en})],20);
              
              % replace with normal
              while ~isempty(rli_removed)
                rl{samples(i)}=replace_rl(rl{samples(i)},rli_removed(1),0);
                %update indices of segments to remove
                rli_removed=find_rl(rl{samples(i)},peeloff_type.regions_rl(regions_rli,1:2)-chr_zero,0); 
              end
              x(pos,samples(i))=0;        
              y(pos,samples(i))=0;
              C2.dat(in_chr(pos),samples(i))=0;
              
              sc(pos)=sc(pos)-regs{k}(n).segments{samples(i)}.vec./ns;
              
              if show_plots>0
                subplot(2,1,1);
                plot(sc); hold on;
                subplot(2,1,2);
                tmp=zeros(length(sc),1);
                tmp(pos)=regs{k}(n).segments{samples(i)}.vec./ns;
                if show_plots==1
                  plot(tmp,'r-');
                  axis(ax);
                  disp([ 'min=' num2str(min(sc)) '; mean=' num2str(mean(regs{k}(n).segments{samples(i)}.vec./ns))]);
                  pause;
                end
              end
              if show_plots==2
                X=reorder_D_rows(C2,in_chr);
                X.dat=x*[1-2*(k==2)];
                X.fxa{1}=sc;
                X.fxa{2}=sc;
                display_D(X,[],[],{'snp_fxa',struct('items',{{{2,5,'fxascore',1,[0 max(orig_sc)],score_thresh(k)}}})});
              end
            end
          end
         
         [mx,mi]=max(sc);
         verbose(['max(score)=' num2str(mx) '(?> ' num2str(score_thresh(k)) ')' ],20);
         if show_plots==1
           figure(1); clf;
           subplot(2,1,1);
           plot(sc); hold on;
           ax=axis;
           pause
         end
       end
       
     end
   end
   
 case 'iterate'
  orig_ads=ads;
  orig_q=q;
  % find most significant peak in amp or del
  for k=1:2
    [mx(k),mxi(k)]=max(-log10(q{k}));
  end
  [mxb,mxbk]=max(mx);
  ch=C.chrn(mxi(mxbk));
  cur_mx_q=q{mxbk}(mxi(mxbk));
  k=mxbk;
  while cur_mx_q<=qv_thresh 
    verbose(['Chr ' num2chromosome(ch)],20);
    
    in_chr=find(C.chrn==ch);
    chr_zero=min(in_chr)-1;
    chr_max = max(in_chr); % Marks last snp on current chromosome
    
    switch score_type.method
        case 'nxa'
            x=C.dat(in_chr,:);
            if k==1 % amp
                th=score_type.amp_thresh;
            else % del
                x=-x;
                th=-score_type.del_thresh;
            end
            x(x<th)=0;
        otherwise
            error('No such score method');
    end
    
    y=(x>th);
    sc=nanmean(x,2);
    orig_sc=sc;
    rl=runlength(y);
    ns=size(y,2);
       
    [mx,mi]=max(sc);
    verbose(['max(score)=' num2str(mx) '(?> ' num2str(score_thresh(k)) ')' ],20);
    if (show_plots==1)
      figure(1); clf;
      subplot(2,1,1);
      plot(sc); hold on;
      ax=axis;
    end
       
    [sads,sadsi]=sort(ads{k});
    sq=log(q{k}(sadsi));
       
    if show_plots==2
      X=reorder_D_rows(C2,in_chr);
      X.dat=x*[1-2*(k==2)];
      X.fxa{1}=sc;
      X.fxa{2}=sc;
      display_D(X,[],[],{'snp_fxa',struct('items',{{{2,5,'fxascore',1,[0 max(orig_sc)],score_thresh(k)}}})});
    end
       
    n=length(regs{k})+1;
         
    % find samples 
    % HANDLE PEAKS WITH MORE THAN ONE SEGMENT!!!
    rg=find(sc==mx);
    if any(diff(rg)>1)
      warning('the peak has more than one segment. should fix enlarged peak');
    end
    samples=find(y(mi,:));
         
    %      CC.sdesc=C.sdesc; CC.dat=C2.dat;
    %      save CC.x.mat CC
    
    % find segments in samples
    rli=find_rl(rl(samples),mi);
    
    regs{k}(n).chrn=ch;
    v=zeros(1,ns);
    v(samples)=1;
    regs{k}(n).samples=v;
    regs{k}(n).peak_st=min(rg)+chr_zero;
    regs{k}(n).peak_en=max(rg)+chr_zero;
    regs{k}(n).peak=round(0.5*(regs{k}(n).peak_st+regs{k}(n).peak_en));
    regs{k}(n).qv=orig_q{k}(in_chr(mi));
    regs{k}(n).resid_qv=exp(interp_pwl(sads,sq,mx));
    regs{k}(n).st=max([0; find(orig_sc(1:mi)<score_thresh(k)) ])+1+chr_zero;
    regs{k}(n).en=min([ find(orig_sc(mi:end)<score_thresh(k)); length(sc)-mi+1 ])-1+mi+chr_zero-1; 
    %% added -1 at end of prev. line since we want the last one that is still above score
    regs{k}(n).score=orig_ads{k}(in_chr(mi));
    verbose(['peak at ' genomic_location(C,{regs{k}(n).peak_st:regs{k}(n).peak_en})],20);
         
    z=zeros(size(x));
    for i=1:length(samples)
      for j=1:length(rli{i})
        rng=rl{samples(i)}(rli{i},1):rl{samples(i)}(rli{i},2);
        z(rng,samples(i))=x(rng,samples(i));
      end
    end
    %      if (k==2) && (ch==14)
    %        keyboard
    %      end
    
    % commented since we want to use the original wide region
    % use next line if you want to use all the data to look for the next boundary
    switch wide_type.method
     case 'robust'
      regs{k}(n).robust_reg_st = max([1; find(orig_sc(1:mi)<(mx- ...
                                                        score_thresh(k)))])+chr_zero;
      right_half = (mi:length(sc))+1;
      regs{k}(n).robust_reg_en = min([ right_half(find(orig_sc(mi:end)< ...
                                                       (mx-score_thresh(k))))+chr_zero-1 max(in_chr)]); 
      [wide_start wide_end thresh range_dist{k}] = robust_peak(regs{k}(n),ads{k},d{k},perm_ads{k},mx-score_thresh(k),score_type,regs{k}(n).robust_reg_st,regs{k}(n).robust_reg_en,conf_level(k),range_dist{k});
      
      regs{k}(n).peak_wide_st = max(chr_zero+1,wide_start-1);
      regs{k}(n).peak_wide_en = min(chr_max,wide_end+1);
     case 'origloo'
      [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en]=add_wide_peak(z,chr_zero);
      %         [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en]=add_wide_peak(regs{k}(n),x,chr_zero);
     case 'leave-k-out'
      %         if (k==2)
      %           keyboard
      %         end
      [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en]=add_wide_peak_k(z,chr_zero,wide_type.k);
     otherwise 
      error('no such method');
    end
    [regs{k}(n).peak_wide_st regs{k}(n).peak_st regs{k}(n).peak_en regs{k}(n).peak_wide_en]
    regs{k}(n).segments=cell(1,ns);
    for i=1:length(samples)
      regs{k}(n).segments{samples(i)}.rli=rli{i};
      regs{k}(n).segments{samples(i)}.st=rl{samples(i)}(rli{i},1)+chr_zero;
      regs{k}(n).segments{samples(i)}.en=rl{samples(i)}(rli{i},2)+chr_zero;
      pos=rl{samples(i)}(rli{i},1):rl{samples(i)}(rli{i},2);
      regs{k}(n).segments{samples(i)}.vec=x(pos,samples(i));
      
      verbose(['removing segment from ' C.sdesc{samples(i)} ' at ' ...
               genomic_location(C,{regs{k}(n).segments{samples(i)}.st:regs{k}(n).segments{samples(i)}.en})],20);
      
      % replace with normal
      rl{samples(i)}=replace_rl(rl{samples(i)},rli{i},Inf); 
      x(pos,samples(i))=Inf;        
      y(pos,samples(i))=Inf;
      C2.dat(in_chr(pos),samples(i))=Inf;
      
      sc(pos)=sc(pos)-regs{k}(n).segments{samples(i)}.vec./ns;
      
      if show_plots>0
        subplot(2,1,1);
        plot(sc); hold on;
        subplot(2,1,2);
        tmp=zeros(length(sc),1);
        tmp(pos)=regs{k}(n).segments{samples(i)}.vec./ns;
        if show_plots==1
          plot(tmp,'r-');
          axis(ax);
          disp([ 'min=' num2str(min(sc)) '; mean=' num2str(mean(regs{k}(n).segments{samples(i)}.vec./ns))]);
          pause;
        end
      end
      if show_plots==2
        X=reorder_D_rows(C2,in_chr);
        X.dat=x*[1-2*(k==2)];
        X.fxa{1}=sc;
        X.fxa{2}=sc;
        display_D(X,[],[],{'snp_fxa',struct('items',{{{2,5,'fxascore',1,[0 max(orig_sc)],score_thresh(k)}}})});
      end
    end
  end
    
    %%%% STOPPED HERE
%     [mx,mi]=max(sc);
%          verbose(['max(score)=' num2str(mx) '(?> ' num2str(score_thresh(k)) ')' ],20);
%          if show_plots==1
%            figure(1); clf;
%            subplot(2,1,1);
%            plot(sc); hold on;
%            ax=axis;
%            pause
%          end
%        end
  
  
 otherwise
  error('no such iter_type');
end
