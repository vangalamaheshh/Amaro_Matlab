function count = count_ziggs_by_sample(Q,thresh)
  
  samples = unique(Q(:,5));
  count = zeros(1,length(samples));
  
  for i=1:length(samples)
    count(i) = length(find(Q(find(Q(:,5)==i),4)>thresh));
  end
  
    
