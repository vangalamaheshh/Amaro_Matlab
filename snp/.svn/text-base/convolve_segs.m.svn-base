function htot = convolve_segs(segs,do_disp,max_bin)
  
  htot = segs{1};
  
  if do_disp
    tic
  end
    
  for j=2:size(segs,2)
    if do_disp
      if mod(j,100) == 0
        disp(j)
      end
    end
    htot = conv(htot,segs{j});
    if length(htot) > max_bin
      htot(max_bin) = sum(htot(max_bin:end));
      htot = htot(1:max_bin);
    end
  end
  
  if do_disp
  toc
  end
  
