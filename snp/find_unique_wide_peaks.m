function unique_wide_peaks = find_unique_wide_peaks(wide_peaks)
  
  size_wide_peaks = wide_peaks(:,2)-wide_peaks(:,1)+1;
  [sx si] =sort(size_wide_peaks);
  unique_wide_peaks = wide_peaks(si,:); %% wide peaks now sorted from smallest
                                 %to largest
  
  for j=1:size(unique_wide_peaks,1)-1
    k=j+1;
    while k<=size(unique_wide_peaks,1)
      %if ~isempty(intersect(unique_wide_peaks(j,1):unique_wide_peaks(j,2),unique_wide_peaks(k,1):unique_wide_peaks(k,2)))
      if unique_wide_peaks(j,1) >= unique_wide_peaks(k,1) & unique_wide_peaks(j,2) <= unique_wide_peaks(k,2) ...
        %% if smaller wide peak contained entirely within larger wide
        %peak
        unique_wide_peaks = unique_wide_peaks(setdiff(1:size(unique_wide_peaks,1),k),:);
      else
        k=k+1;
      end
    end
    
  end    
      
    
  
