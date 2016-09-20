function htot = conv_n_times(hc,n)
  
  %% Given a vector hc, conv_n_times computes the convolution of hc with
  %itself a total of n times in a relatively efficient manner
    
    temp_h = cell(1,n);
    temp_h = cell(1,n);
    temp_h{1} = hc;
    
    bin_string = dec2bin(n);
    
    if str2num(bin_string(end)) == 1
      htot = hc;
    else
      htot = 1;
    end
    
    for j=1:length(bin_string)-1 %% j represents the convolution of
                                 %2^(j-1) copies of hc with itself
      disp(['Round ' num2str(j) ' of ' num2str(length(bin_string)-1)]);
      cur_h = temp_h{2^(j-1)};
      temp_h{2^j} = conv(cur_h,cur_h);
      

      if str2num(bin_string(end-j)) == 1
        disp('Adding result to total convolution!');
        htot = conv(htot,temp_h{2^j});
      end
    end
    
        
    
