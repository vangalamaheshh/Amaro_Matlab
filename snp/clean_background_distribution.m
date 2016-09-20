function [regs new_d] = clean_background_distribution(D,Z,QA,QD,ads,d,p,score_type,score_thresh,clean_thresh,plot_diff)
  
  if ~exist('plot_diff','var') || isempty(plot_diff)
    plot_diff = 1;
  end
  
  disp('Performing background clean-up!')
  
  bonf_n = size(Z.dat{1},1);
  regs = cell(1,2);
  cleared = cell(1,2);
  Qamp = QA;
  Qdel = QD;
  
  for k=1:2
    cleared{k} = zeros(size(Z.dat{1},1),size(Z.dat{1},2),'single');
    t{k}=flipud(cumsum(flipud(d{k})));
    if k==1
      cleared{k}(D.dat < -0.1) = 1;
    else
      cleared{k}(D.dat > 0.1) = 1;
    end
    [min_p mi] = min(p{k});
    if k==1
      Q = QA;
    else
      Q = QD;
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
      n=length(regs{k})+1; % n = current reg number
      rg=find(p{k}(in_chr)==min_p);
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
      regs{k}(n).bonfp = bonf_p; % bonferonni p-value of alteration
      regs{k}(n).st=max([0; find(sc(1:mi)<score_thresh(k)) ...
                        ])+1+chr_zero; 
      regs{k}(n).en=min([ find(sc(mi:end)<score_thresh(k)); length(sc)-mi+1 ...
                        ])-1+mi+chr_zero-1; % finds last snp on chr
                                            % scoring above thresh
     regs{k}(n).score=mx; % score of peak
       
     if k==1
       disp(['Cleaning amplification peak at ' genomic_location(Z,{regs{k}(n).peak_st:regs{k}(n).peak_en})]);
     else
       disp(['Cleaning deletion peak at ' genomic_location(Z,{regs{k}(n).peak_st:regs{k}(n).peak_en})]);
     end
     
     % Now peel-off & clean-up
     
     [Znew,Qnew,rm_segments] = ziggurat_peeloff(Z.dat{k},Q,regs{k}(n).peak);
     
     Qrm = Q(rm_segments,:);
          
     rm_samples = unique(Qrm(:,5));
     nonrm_samples = setdiff(1:ns,rm_samples); 
          
     if k==1
       Q = Qnew;
     else
       Q = Qnew;
     end
     
     rm_length = zeros(1,length(rm_samples));
     
     for i=1:length(rm_samples)
       cur_segs = find(Qrm(:,5) == rm_samples(i));
       rm_start = min(Qrm(cur_segs,2));
       rm_end = max(Qrm(cur_segs,3));
       rm_length(i) = rm_end-rm_start+1;
       cleared{k}(rm_start:rm_end,rm_samples(i)) = 1;
     end
     
     mean_length = mean(rm_length);
     mean_start = max(1,round(regs{k}(n).peak-mean_length/2));
     mean_end = min(size(Znew,1),round(regs{k}(n).peak+mean_length/2));
     
     cleared{k}(mean_start:mean_end,nonrm_samples) = 1;
          
     Z.dat{k} = Znew;
     
     %% Recompute ads and p-values
     
     ads{k} = nanmean(Z.dat{k},2);
     p{k} = t{k}(min(1+floor(ads{k}/score_type.res),length(t{k})));
          
     % Find new bonf_p,max score and chromosome 
     [min_p mi] = min(p{k});
     bonf_p = min_p*bonf_n
     mx = ads{k}(mi);
     ch = Z.chrn(mi);
     
    end
  end
  
  disp('Recomputing cleaned background distributions!');
  
  cleared_snps = cell(1,2);
  
  for k=1:2
    cleared_snps{k} = sum(cleared{k});
  end
  
  
  focalsA = zeros(size(Z.dat{1},1),size(Z.dat{1},2));
  focalsD = focalsA;
  
  for i=1:size(Qamp,1)
    focalsA(Qamp(i,2):Qamp(i,3),Qamp(i,5)) = focalsA(Qamp(i,2):Qamp(i,3), ...
                                                     Qamp(i,5))+Qamp(i,4);
  end
  
  for i=1:size(Qdel,1)
    focalsD(Qdel(i,2):Qdel(i,3),Qdel(i,5)) = focalsD(Qdel(i,2):Qdel(i,3), ...
                                                     Qdel(i,5))+Qdel(i,4);
  end
  
  
  [new_q{1},new_p{1},new_d{1},ads{1}] = zigg_score_permutations(focalsA, ...
                                                    score_type.res,cleared_snps{1});

  [new_q{2},new_p{2},new_d{2},ads{2}] = zigg_score_permutations(focalsD, ...
                                                    score_type.res,cleared_snps{2});
  
  if plot_diff
    figure(1)
    title('Amplification d')
    plot(d{1},'Color','blue');
    hold on
    plot(new_d{1},'Color','red');
    legend('Old d','New d');
    
    figure(2)
    title('Deletion d')
    plot(d{2},'Color','blue');
    hold on
    plot(new_d{2},'Color','red');
    legend('Old d','New d');
    
  end
