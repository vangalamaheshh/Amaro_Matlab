function [local_max max_vals] = find_local_score_max(scores,range)
  
  local_max = [];
  
  range_start = range(1);
  range_end = range(2);
  
  data_range = scores(range_start:range_end);
  
  bpts = find(diff(data_range) ~= 0);
  bpts = union(bpts,range_end-range_start);
  
  diff_bpts = diff(data_range(bpts));
  
  if diff_bpts(1) < 0 
    local_max = cat(1,local_max,[0 bpts(1)]);
  end
  
  for j=1:length(diff_bpts)
    if (diff_bpts(j) > 0 & (j == length(diff_bpts) | diff_bpts(j+1) < 0))
      if j < length(diff_bpts)
        local_max = cat(1,local_max,[bpts(j)+1 bpts(j+1)
]);
      else
        local_max = cat(1,local_max,[bpts(j)+1 bpts(j+1)]);
      end
    end
  end
  
  max_vals = data_range(round(mean(local_max,2)));
  [max_vals si] = sort(max_vals,'descend');
  local_max = local_max(si,:);
  
  local_max(:,1) = local_max(:,1)+range_start;
  local_max(:,2) = local_max(:,2)+range_start-1;
  
  
