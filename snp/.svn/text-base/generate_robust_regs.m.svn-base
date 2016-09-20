function regs = generate_robust_regs(Z,QA,QD,ads,nperm_ads,d,q,score_type,score_thresh,conf_level)
  
  regs = cell(1,2);
      
  for k=1:2
    if k==1
      Q=QA;
    else
      Q=QD;
    end
    
    orig_ads = ads{k};
    perm_ads = nperm_ads{k};
    bg_dist = d{k};
    bpps = length(Q)/size(Z.dat{k},1);
      
    chr_lengths = zeros(1,max(Z.chrn));
    for ch=1:max(Z.chrn)
      chr_lengths(ch) = length(find(Z.chrn==ch));
    end
    max_chr_snps = max(chr_lengths);
    
    range_dist = cell(1,max_chr_snps);
    
    for ch=1:max(Z.chrn)
      in_chr=find(Z.chrn==ch); % Find snps within that chromosome
      chr_zero= min(in_chr)-1; % Marks first snp not on current chromosome
      chr_max = max(in_chr); % Marks last snp on current chromosome
      x=Z.dat{k}(in_chr,:); % x contains only data from present
                            % chromosome
      y=(x>0);
      sc=nanmean(x,2); % take row-wise mean of each data
      orig_sc = sc;
      ns=size(x,2); % ns = # samples
      [mx mi] = max(sc);
      
      while mx >= score_thresh(k)
        rg=find(sc==mx);
        
        if any(diff(rg)>1)
          warning(['the peak has more than one segment']);
          segs = [0 find(diff(rg)~=1) length(rg)];
          [sx si] = max(diff([0 find(diff(rg)>1) length(rg)]));
          rg = rg(segs(si)+1:segs(si+1));  
          mi=rg(1);
        end
        
        right_half = (mi:length(sc))+1;
        
        robust_reg_start = max([1; find(sc(1:mi)<(mx-score_thresh(k))) ...
                   ])+chr_zero; 
        
        robust_reg_end = min([ right_half(find(sc(mi:end)< (mx-score_thresh(k))))+chr_zero-1 max(in_chr)]);  
        
        
        % See if there is more than one local maximum in robust region
        %local_max = find_local_score_max(ads{k},[robust_reg_start ...
        %                    robust_reg_end]);
        
        local_max = [min(rg)+chr_zero max(rg)+chr_zero];
        
        unique_wide_peaks = [];
        if size(local_max,1) > 1
          wide_peaks = [];
          for i=1:size(local_max,1)
            disp(['Testing local max ' num2str(i) ' of ' num2str(size(local_max,1))]);
            temp_reg.score = sc(local_max(i,1)-chr_zero);
            temp_reg.peak_st = local_max(i,1);
            temp_reg.peak_en = local_max(i,2);
            [wide_peaks(i,1) wide_peaks(i,2) thresh1] = robust_peak(temp_reg,ads{k},bg_dist,perm_ads,conf_level(k),mx-score_thresh(k),score_type,robust_reg_start,robust_reg_end,bpps,range_dist,1);
            
            %[wide_peaks(i,1) wide_peaks(i,2) thresh2 range_dist] = robust_peak(temp_reg, ...
%                                                  ads{k},bg_dist,perm_ads,conf_level(k),thresh1,score_type,wide_start1,wide_end1,bpps,range_dist,2);
            %[wide_peaks(i,1) wide_peaks(i,2)] =robust_on_segments(Q, ...
            %                                                  temp_reg, ...
            %                                                  ads{k},amp_hist{k},res,score_thresh(k),wide_start1,wide_end1,conf_level(k),nperm);
            %if i>1 & wide_peaks(i,1) <= wide_peaks(i-1,1) & wide_peaks(i,2) ...
            %      >= wide_peaks(i-1,2)
            %  disp('Local max are not distinct.  Stopping search.');
            %  break;
            %end
            
          end
          unique_wide_peaks = find_unique_wide_peaks(wide_peaks);
        else
          %disp(['Testing local max ' num2str(1) ' of ' num2str(size(local_max,1))]);
          temp_reg.score = sc(rg(1));
          temp_reg.peak_st = min(rg)+chr_zero;
          temp_reg.peak_en = max(rg)+chr_zero;
          [unique_wide_peaks(1,1) unique_wide_peaks(1,2) thresh1 range_dist] = robust_peak(temp_reg,ads{k}, ...
                                                        bg_dist,perm_ads,conf_level(k),mx-score_thresh(k),score_type,robust_reg_start,robust_reg_end,bpps,range_dist,2);
 %         [unique_wide_peaks(1,1) unique_wide_peaks(1,2) thresh2 range_dist] = robust_peak(temp_reg,ads{k},bg_dist,perm_ads,conf_level(k),thresh1,score_type,wide_start1,wide_end1,bpps,range_dist,2);
          %[unique_wide_peaks(1,1) unique_wide_peaks(1,2)] = robust_on_segments(Q,temp_reg,ads{k},amp_hist{k},res,score_thresh(k),wide_start1,wide_end1,conf_level(k),nperm);
        end
        
        if size(unique_wide_peaks,1) > 1
          disp(['Found more than one independent wide peak in robust ' ...
                'region!']);
        end
        
        for i=1:size(unique_wide_peaks) %% for each unique wide
                                        %peak, generate a reg
          n=length(regs{k})+1; % n = current reg number
          [mx mi] = max(sc(unique_wide_peaks(i,1)-chr_zero: ...
                           unique_wide_peaks(i,2)-chr_zero));
          rg = find(sc == mx);
          mi = mi+unique_wide_peaks(i,1)-chr_zero;
          samples = find(y(mi,:));
          regs{k}(n).chrn=ch; % sets chr number of current reg
          v=zeros(1,ns); 
          v(samples)=1; % sets samples with alteration to 1
          regs{k}(n).samples=v; 
          regs{k}(n).peak_st = min(rg)+chr_zero;
          regs{k}(n).peak_en = max(rg)+chr_zero;
          regs{k}(n).peak=round(0.5*(regs{k}(n).peak_st+regs{k}(n).peak_en));
          regs{k}(n).qv = q{k}(in_chr(mi)); 
          right_half = (mi:length(sc))+1;
          regs{k}(n).st=max([0; find(orig_sc(1:mi)<score_thresh(k)) ...
                            ])+chr_zero+1; % finds first snp on chr
                                           % scoring above score_thresh
          regs{k}(n).en=min([ right_half(find(orig_sc(mi:end)< ...
                                              score_thresh(k)))+ ...
                              chr_zero+1 max(in_chr)]); % finds last snp on chr scoring above thresh  
          regs{k}(n).score = mx;
          
          %[regs{k}(n).peak_wide_st regs{k}(n).peak_wide_en] = robust_peak(temp_reg,orig_ads,bg_dist,perm_ads,conf_level(k),mx-score_thresh(k),score_type,unique_wide_peaks(i,1),unique_wide_peaks(i,2),bpps,range_dist,2);
          %[regs{k}(n).peak_wide_st regs{k}(n).peak_wide_en] = robust_on_segments(Q,regs{k}(n),ads{k},d{k},size(Z.dat{k},2),score_type,robust_reg_start,robust_reg_end,conf_level(k),nperm,num_iters);
          %[regs{k}(n).peak_wide_st regs{k}(n).peak_wide_en] = robust_on_segments(Q,regs{k}(n),ads{k},d{k},size(Z.dat{k},2),score_type,unique_wide_peaks(i,1),unique_wide_peaks(i,2),conf_level(k),nperm,num_iters);
          %regs{k}(n).peak_wide_st = max(chr_zero+1, ...
        %                                regs{k}(n).peak_wide_st-1);
          %regs{k}(n).peak_wide_en = min(chr_max,regs{k}(n).peak_wide_en+1);
          regs{k}(n).peak_wide_st = max(chr_zero+1,unique_wide_peaks(i,1)-1);
          regs{k}(n).peak_wide_en = min(chr_max,unique_wide_peaks(i,2)+1);
          
          if k==1
            disp(['Amplification wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
          else
            disp(['Deletion wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
          end
          %% peel-off segments
          [Z.dat{k} Q] = ziggurat_peeloff(Z.dat{k},Q,regs{k}(n).peak);
        end
        
        %% recalculate score  
        ads{k} = nanmean(Z.dat{k},2);
        x = Z.dat{k}(in_chr,:);
        y = (x>0);
        sc = nanmean(x,2);
        [mx mi] = max(sc);
        
      end %% end of while loop
      
    end %% end of ch loop
    
  end %% end of k=1:2 loop
  
