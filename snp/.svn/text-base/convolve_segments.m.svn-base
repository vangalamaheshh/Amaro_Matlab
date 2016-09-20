function htot = convolve_segments(seg_hists,do_disp,max_bin)
  
  if ~exist('do_disp') || isempty(do_disp)
    do_disp = 0;
  end
  
  ha_cur = seg_hists;
    
  for i=1:ceil(log2(size(seg_hists,2)))
    if do_disp
      disp(['Round ' num2str(i) ' of ' num2str(ceil(log2(size(seg_hists, ...
                                                        2))))]);
      tic
    end
    ha_new = cell(1,ceil(size(ha_cur,2)/2));
    for j=1:floor(size(ha_cur,2)/2)
      if do_disp
        if mod(j,1000) == 0
          disp(j)
        end
      end
      ha_new{j} = sparse_conv(ha_cur{2*j-1},ha_cur{2*j});
    end
    if size(ha_cur,2) > 2*j
      ha_new{j+1} = ha_cur{2*j+1};
    end
    if do_disp
       toc;
    end
    ha_cur = ha_new;
  end
  
  htot = ha_cur{1};
  
  if exist('max_bin') && ~isempty(max_bin)
    if length(htot) > max_bin
      htot(max_bin) = sum(htot(max_bin:end));
      htot = htot(1:max_bin);
    end
  end
