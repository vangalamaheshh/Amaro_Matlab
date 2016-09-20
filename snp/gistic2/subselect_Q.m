function Qs = subselect_Q(Qs,sample_idx)
    
  fields = intersect({'amp','del','aod','doa'},fieldnames(Qs));
    
  new_Q = struct();
  
  for i=1:length(sample_idx)
    sample_lookup(sample_idx(i)) = i;
  end


  for k=1:length(fields)
    cur_Q = getfield(Qs,fields{k});
    cur_idx = [];
    for i=1:length(sample_idx)
      if mod(i,100) == 0
        disp([k i])
      end
      new_idx = find(cur_Q(:,5) == sample_idx(i));
      cur_idx = [cur_idx; new_idx];
    end
    temp_Q = cur_Q(cur_idx,:);
    temp_Q(:,5) = sample_lookup(temp_Q(:,5));
    new_Q = setfield(new_Q,fields{k},temp_Q);
  end
  
  new_Q.sdesc = Qs.sdesc(sample_idx,:);
  new_Q.header = Qs.header;
  
  Qs = new_Q;
  
  
  
  
  
