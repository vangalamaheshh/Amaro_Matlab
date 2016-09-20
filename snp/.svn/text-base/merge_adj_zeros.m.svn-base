function M = merge_adj_segs(M)
  
  % Takes a matrix M corresponding to chromosomal segments and merges
  % adjacent rows whose value in column 4 are equal
    
    i=1;
    while i < size(M,1)
      if M(i,4) == M(i+1,4)
        M(i,3) = M(i+1,3);
        M=M(setdiff(1:size(M,1),i+1),:);
      else
        i=i+1;
      end
    end
    
        
