function [regs new_d] = generate_regs_by_zigg_peel_off(Z,QA,QD,ads,d,p,score_type,score_thresh,conf_level,two_stage,iter_cleanup,clean_thresh)
  
  if ~exist('two_stage','var') || isempty(two_stage)
    two_stage = 0;
  end
      
  if ~exist('iter_cleanup','var') || isempty(iter_cleanup)
    iter_cleanup = 0;
  end
  
  if ~exist('clean_thresh','var') || isempty(clean_thresh)
    clean_thresh = 0.05;
  end
  
    if ~iter_cleanup
      merged_p = cat(1,p{1},p{2});
      merged_q = calc_fdr_value(merged_p);
      q{1} = merged_q(1:length(p{1}));
      q{2} = merged_q(length(p{1})+1:end);
      
      disp('Not using iterative clean-up!')
      regs = cell(1,2);
      for k=1:2
        if k==1
          Q=QA;
        else
          Q=QD;
        end
      for ch = 1:max(Z.chrn)
        in_chr=find(Z.chrn==ch); % Find snps within that chromosome
        chr_zero=min(in_chr)-1; % Marks first snp not on current chromosome
        chr_max = max(in_chr); % Marks last snp on current chromosome
        x=Z.dat{k}(in_chr,:); % x contains only data from present chromosome
        y=(x>0);
        sc=nanmean(x,2); % take row-wise mean of each data
        orig_sc = sc;
        ns=size(x,2); % ns = # samples
        [mx,mi]=max(sc); % mx = maximum score on chr, mi = index of first maximum
                         % scoring snp
        
        [sads,sadsi]=sort(ads{k}); %sads contains g-scores in sorted order
                                   %(lowest to highest), sadsi is index
                                   %permutation corresponding to this ordering
        sq=log(q{k}(sadsi)); %sq is log of q-values in sorted order
        
        while mx >= score_thresh(k)
          n=length(regs{k})+1; % n = current reg number
          rg=find(sc==mx);
          if any(diff(rg)>1)
            warning(['the peak has more than one segment']);
            segs = [0 find(diff(rg)~=1)' length(rg)];
            [sx si] = max(diff([0 find(diff(rg)>1)' length(rg)]));
            rg = rg(segs(si)+1:segs(si+1));  
            mi=rg(1);
          end
          samples = find(y(mi,:));
          regs{k}(n).chrn=ch; % sets chr number of current reg
          v=zeros(1,ns); 
          v(samples)=1; % sets samples with alteration to 1
          regs{k}(n).samples=v; 
          regs{k}(n).peak_st=min(rg)+chr_zero; % first snp in max. segment
          regs{k}(n).peak_en=max(rg)+chr_zero; % last snp in max. segment
          regs{k}(n).peak=round(0.5*(regs{k}(n).peak_st+regs{k}(n).peak_en));% middle snp of peak region
          regs{k}(n).qv=q{k}(in_chr(mi)); % q-value of alteration
          regs{k}(n).resid_qv=exp(interp_pwl(sads,sq,mx)); % resid-qv of alteration
          regs{k}(n).st=max([0; find(orig_sc(1:mi)<score_thresh(k)) ...
                            ])+1+chr_zero; 
          % finds first snp on chr scoring above score_thresh
          regs{k}(n).en=min([ find(orig_sc(mi:end)<score_thresh(k)); length(sc)-mi+1 ...
                            ])-1+mi+chr_zero-1; % finds last snp on chr
                                                % scoring above thresh
          robust_reg_start = max([0; find(sc(1:mi)<(mx-score_thresh(k))) ...
                            ])+1+chr_zero; 
          robust_reg_end = min([ find(sc(mi:end)<(mx-score_thresh(k))); length(sc)-mi+1 ...
                            ])-1+mi+chr_zero-1; % finds last snp on chr
                                                % scoring above thresh
          regs{k}(n).score=mx; % score of peak
                  
          bpps = length(Q)/size(Z.dat{k},1);
          if ~two_stage
            [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en] = ...
                add_wide_peak_conf(regs{k}(n),chr_zero,chr_max,ads{k}, ...
                                   d{k},conf_level(k),score_thresh(k),score_type,bpps,0);
            if k==1
              disp(['Amplification wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
            else
              disp(['Deletion wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
            end
          else
            new_d = xcorr(d{k},d{k});
            new_d = new_d(median(1:length(new_d)):end);
            [regs{k}(n).wide_peaks peak_wide_st peak_wide_en] = robust_peak(regs{k}(n),ads{k},new_d, ...
                                                conf_level(k), ...
                                                              mx-score_thresh(k),score_type,bpps,robust_reg_start,robust_reg_end);
           regs{k}(n).peak_wide_st = max(chr_zero+1,peak_wide_st-1);
           regs{k}(n).peak_wide_en = min(chr_max,peak_wide_en+1);
            if k==1
              disp(['Amplification wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
            else
              disp(['Deletion wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
            end
          end
                              
          % Now peel-off
          
          [Z.dat{k},Q] = ziggurat_peeloff(Z.dat{k},Q,regs{k}(n).peak);
          x = Z.dat{k}(in_chr,:);
          y = (x>0);
          sc = nanmean(x,2);
          [mx,mi]=max(sc);
          
        end
      end
    end
    new_d = d;
  else %% Using iterative cleanup
    disp('Using iterative clean-up!')
    bonf_n = 2*size(Z.dat{1},1);
    %ampD = [d{1}];
    %delD = [d{2}];
    regs = cell(1,2);
    clearedA = zeros(size(Z.dat{1},1),size(Z.dat{1},2),'single');
    clearedD = clearedA;
    % Find min p-value (on either amp or del)
    [minp1 mi1] = min(p{1});
    [minp2 mi2] = min(p{2});
    if minp1 <= minp2
      min_p = minp1;
      mi = mi1;
      k=1;
      Q=QA;
    else
      min_p = minp2;
      mi = mi2;
      k=2;
      Q=QD;
    end
    % Find max score and chromosome
    bonf_p = min_p*bonf_n
    mx = ads{k}(mi);
    ch = Z.chrn(mi);
    while bonf_p <= clean_thresh
      in_chr=find(Z.chrn==ch); % Find snps within that chromosome
      chr_zero=min(in_chr)-1; % Marks first snp not on current chromosome
      chr_max = max(in_chr); % Marks last snp on current chromosome
      mi = mi-chr_zero;
      x=Z.dat{k}(in_chr,:); % x contains only data from present chromosome
      y=(x>0);
      sc=nanmean(x,2); % take row-wise mean of each data
      ns=size(x,2); % ns = # samples
      [sads,sadsi]=sort(ads{k}); %sads contains g-scores in sorted order
                                   %(lowest to highest), sadsi is index
                                   %permutation corresponding to this ordering
      % FIXME sq=log(q{k}(sadsi)); %sq is log of q-values in sorted order
      n=length(regs{k})+1; % n = current reg number
      rg=find(sc==mx);
      if any(diff(rg)>1)
        warning(['the peak has more than one segment']);
        segs = [0 find(diff(rg)~=1) length(rg)];
        [sx si] = max(diff([0 find(diff(rg)>1) length(rg)]));
        rg = rg(segs(si)+1:segs(si+1));  
        mi=rg(1);
      end
      samples = find(y(mi,:));
      regs{k}(n).chrn=ch; % sets chr number of current reg
      v=zeros(1,ns); 
      v(samples)=1; % sets samples with alteration to 1
      regs{k}(n).samples=v; 
      regs{k}(n).peak_st=min(rg)+chr_zero; % first snp in max. segment
      regs{k}(n).peak_en=max(rg)+chr_zero; % last snp in max. segment
      regs{k}(n).peak=round(0.5*(regs{k}(n).peak_st+regs{k}(n).peak_en));% middle snp of peak region
      regs{k}(n).bonfp = bonf_p; % q-value of alteration
      % FIXME!!! regs{k}(n).resid_qv=exp(interp_pwl(sads,sq,mx)); % resid-qv of alteration
      regs{k}(n).st=max([0; find(sc(1:mi)<score_thresh(k)) ...
                        ])+1+chr_zero; 
      % finds first snp on chr scoring above score_thresh
      regs{k}(n).en=min([ find(sc(mi:end)<score_thresh(k)); length(sc)-mi+1 ...
                        ])-1+mi+chr_zero-1; % finds last snp on chr
                                            % scoring above thresh
      regs{k}(n).score=mx; % score of peak
      bpps = length(Q)/size(Z.dat{k},1);
      if ~two_stage
        [regs{k}(n).peak_wide_st regs{k}(n).peak_wide_en] = add_wide_peak_conf(regs{k}(n),chr_zero,chr_max,ads{k},d{k},conf_level(k),score_thresh(k),score_type,bpps,0);
      else
        % Stage 1:
        [step1_start,step1_end] = ...
            add_wide_peak_conf(regs{k}(n),chr_zero,chr_max,ads{k},d{k},conf_level(k),score_thresh(k),score_type,bpps,0);
        % Compute new distribution
        new_d = xcorr(d{k},d{k});
        new_d = new_d(median(1:length(new_d)):end);
        [regs{k}(n).peak_wide_st,regs{k}(n).peak_wide_en] = ...
            add_wide_peak_conf(regs{k}(n),step1_start-1,step1_end,ads{k}, ...
                               new_d,conf_level(k),score_thresh(k),score_type,bpps,0);
      end
      if k==1
        disp(['Amplification wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
      else
        disp(['Deletion wide peak at ' genomic_location(Z,{regs{k}(n).peak_wide_st:regs{k}(n).peak_wide_en})]);
      end
        
      % Now peel-off & clean-up
      
      rm_segments = find(Q(:,2)<=regs{k}(n).peak_st & Q(:,3)>= ...
                         regs{k}(n).peak_en);
      Qrm = Q(rm_segments,:);
      rm_samples = unique(Qrm(:,5));
      nonrm_samples = setdiff(1:ns,rm_samples); 
      [Znew,Qnew] = ziggurat_peeloff(Z.dat{k},Q,regs{k}(n).peak);
      if k==1
        QA = Qnew;
      else
        QD = Qnew;
      end
      rm_length = zeros(1,length(rm_samples));
      [clear_i clear_j] = find(Z.dat{k}-Znew ~= 0);
      for i=1:length(rm_samples)
        cur_snps = find(clear_j == rm_samples(i));
        rm_start = min(clear_i(cur_snps));
        rm_end = max(clear_i(cur_snps));
        rm_snps = clear_i(cur_snps(find(Znew(rm_start:rm_end,rm_samples(i)) ...
                                        == 0)));
        rm_length(i) = length(rm_snps);
        clearedA(rm_snps,rm_samples(i)) = 1;
        clearedD(rm_snps,rm_samples(i)) = 1;
      end
      mean_length = mean(rm_length); 
      mean_start = max(1,round(regs{k}(n).peak-mean_length/2));
      mean_end = min(size(Znew,1),round(regs{k}(n).peak+mean_length/2));
      if k == 1
        clearedA(mean_start:mean_end,nonrm_samples) = 1;
      else
        clearedD(mean_start:mean_end,nonrm_samples) = 1;
      end
      clearedA_snps = sum(clearedA,1);
      clearedD_snps = sum(clearedD,1);
      Z.dat{k} = Znew;
      
      %% Recompute background distributions for amplifications and
      %deletions!
            
      % Recompute p-values (for new peeled-off data using original d)
      
      for j=1:2
        t{j}=flipud(cumsum(flipud(d{j})));
        ads{j} = nanmean(Z.dat{j},2);
        p{j} = t{j}(min(1+floor(ads{j}/score_type.res),length(t{j})));
      end
      
%%      if k==1
%%        ext_D = cat(1,new_d{1},zeros(size(ampD,1)-length(new_d{1}),1));
%%        ampD = cat(2,ampD,ext_D);
%%      else
%%        ext_D = cat(1,new_d{2},zeros(size(delD,1)-length(new_d{2}),1));
%%        delD = cat(2,delD,ext_D);
%%      end
      
      % Find min p-value (on either amp or del)
       [minp1 mi1] = min(p{1});
       [minp2 mi2] = min(p{2});
       if minp1 <= minp2
         min_p = minp1;
         mi = mi1;
         k=1;
         Q=QA;
       else
         min_p = minp2;
         mi = mi2;
         k=2;
         Q=QD;
       end
       % Find new score_thresh,max score and chromosome 
       bonf_p = min_p*bonf_n
       mx = ads{k}(mi);
       ch = Z.chrn(mi);
       
    end
    disp('Recomputing cleaned background distributions!');

    [new_q{1},new_p{1},new_d{1},ads{1}] = zigg_score_permutations(Z.dat{1},score_type.res,clearedA_snps);
    
    [new_q{2},new_p{2},new_d{2},ads{2}] = zigg_score_permutations(Z.dat{2},score_type.res,clearedD_snps);

    figure(1)
    title('Amplification d')
    plot(log(d{1}),'Color','blue');
    hold on
    plot(log(new_d{1}),'Color','red');
    figure(2)
    title('Deletion d')
    plot(log(d{1}),'Color','blue');
    hold on
    plot(log(new_d{1}),'Color','red');
    
    end
  
