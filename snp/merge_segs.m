function Q = merge_segs(rl,ids)
  
  rlis = rl(ids,1:2);
  Q = [];
  if length(ids) < 1
    return
  else
    diff_ids = diff(ids);
    while length(find(diff_ids == 1)) > 0
      first_id = min(find(diff_ids == 1));
      rlis(first_id,:) = [rlis(first_id,1) rlis(first_id+1,2)];
      rlis = rlis(setdiff(1:size(rlis,1),first_id+1),:);
      ids = setdiff(ids,ids(first_id));
      diff_ids = diff(ids);
    end
    Q = rlis;
  end   
