function amp_hist = gen_amp_histograms(QA,QD,res,by_snps)
  
  if ~exist('by_snps','var') || isempty(by_snps)
    by_snps = 1;
  end

  amp_hist = cell(1,2);
  
  if ~by_snps
  
    amp_hist{1} = histc(QA(:,4),0:res:max(QA(:,4)));
    amp_hist{2} = histc(QD(:,4),0:res:max(QD(:,4)));
  
      
  else
    amp_hist{1} = hist_by_snps(QA(:,4),QA(:,3)-QA(:,2)+1,0:res:max(QA(:, ...
                                                      4)));
    amp_hist{2} = hist_by_snps(QD(:,4),QD(:,3)-QD(:,2)+1,0:res:max(QD(:, ...
                                                      4)));
  end

  amp_hist{1} = amp_hist{1}/sum(amp_hist{1});
  amp_hist{2} = amp_hist{2}/sum(amp_hist{2});
