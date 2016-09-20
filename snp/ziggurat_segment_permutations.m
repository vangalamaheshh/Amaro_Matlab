function [q,p,d,ads] = ziggurat_segment_permutations(M,Q,score_type,seg_assignments,old_d)
  
  res = score_type.res;
      
  if ~exist('seg_assignments','var') || isempty(seg_assignments)
    clean_up = 0;
  else
    disp('Performing clean-up...');
    clean_up = 1;
    gmean = 0;
    for i=1:length(old_d)
      gmean = gmean + i*old_d(i);
    end
    gmean = gmean*res;
  end
    
  n=size(M,2);
  s=size(M,1);   
  
  seg_lengths = Q(:,3)-Q(:,2)+1;
  ads = nanmean(M,2);

  max_amp_bin = ceil(max(max(M))/(n*res));
  ha = cell(1,n);
  disp('Calculating segment histograms for each sample...')
  for i=1:n
    if mod(i,100) == 0
      disp(i)
    end
    cur_segs = find(Q(:,5) == i);
    if isempty(cur_segs)
      ha{i} = [1];
    else
      temp_hists = cell(1,length(cur_segs));
      for j=1:length(cur_segs)
        cur_hist = sparse(1,length(0:res:Q(cur_segs(j),4)/n));
        if ~clean_up
          cur_hist(ceil(Q(cur_segs(j),4)/n/res)) = seg_lengths(cur_segs(j))/s;
          cur_hist(1) = 1-seg_lengths(cur_segs(j))/s;
          temp_hists{j} = cur_hist;
        else
          assigned_regs = find(seg_assignments(cur_segs(j),:));
          if ~isempty(assigned_regs)
            cur_hist(ceil(Q(cur_segs(j),4)/n/res)) = ...
                seg_lengths(cur_segs(j))/s*gmean^(length(assigned_regs))/ ...
                sum(sum(seg_assignments(:,assigned_regs),1));
            cur_hist(1) = 1-cur_hist(ceil(Q(cur_segs(j),4)/n/res));
            temp_hists{j} = cur_hist;
          else
            cur_hist(ceil(Q(cur_segs(j),4)/n/res)) = seg_lengths(cur_segs(j))/s;
            cur_hist(1) = 1-seg_lengths(cur_segs(j))/s;
            temp_hists{j} = cur_hist;
          end
        end
      end
      ha{i} = convolve_segments(temp_hists,0,max_amp_bin);
      %if length(ha{i}) > max_amp_bin
      %  ha{i}(max_amp_bin) = sum(ha{i}(max_amp_bin:end));
      %  ha{i} = ha{i}(1:max_amp_bin);
      %end
    end
  
  end
  
  max_ads_bin = ceil(max(ads)/res+1);
  disp('Computing background distribution!')
  htot = convolve_segments(ha,1,max_ads_bin);
  
  disp('Computing stats...')
  d = htot;
  t = flipud(cumsum(flipud(d)));
  p = t(min(1+floor(ads/res),length(t)));
  q = calc_fdr_value(p);
