function perm_ads = generate_permuted_scores(M,ts,nperm,ampdel)
  
  
  if ~exist('ampdel','var') || isempty(ampdel)
    ampdel = [1 1];
  end
  
  x = runlength(M);
  perm_ads = cell(1,2);
  ns = size(x,2);
  
  
  Qs = cell(1,2);
  
  for j=1:ns
    cur_rl = x{j};
    if ampdel(1)
      amp_ids = find(cur_rl(:,3) > ts(1));
      amp_segs = merge_segs(cur_rl,amp_ids);
      Qs{1} = [Qs{1}; [amp_segs repmat(j,size(amp_segs,1),1)]];
    end
    if ampdel(2)
      del_ids = find(cur_rl(:,3) < -1*ts(2));
      del_segs = merge_segs(cur_rl,del_ids);
      Qs{2} = [Qs{2}; [del_segs repmat(j,size(del_segs,1),1)]];
    end
  end
  
  
  for k=1:2
    if ampdel(k)
      new_ads = cell(1,nperm);
      for j=1:nperm
        disp(['Calculating permutation ' num2str(j) ' of ' num2str(nperm)]);
        shift = floor(size(M,1)*rand(size(Qs{k},1),1));
        temp_Q = Qs{k};
        orig_Q = Qs{k};
        temp_Q(:,1:2) = mod(temp_Q(:,1:2)+repmat(shift,1,2),size(M,1))+1;
        cur_ads = zeros(1,size(M,1));
        for i=1:size(temp_Q,1)
          data = M(orig_Q(i,1):orig_Q(i,2),orig_Q(i,3))';
          if temp_Q(i,1) <= temp_Q(i,2)
            cur_ads(temp_Q(i,1):temp_Q(i,2)) = cur_ads(temp_Q(i,1): ...
                                                       temp_Q(i,2))+data/ns;
          else
            first_half = temp_Q(i,1):size(M,1);
            cur_ads(first_half) = cur_ads(first_half)+data(1:length(first_half))/ns;
            second_half = 1:temp_Q(i,2);
            cur_ads(second_half) = cur_ads(second_half)+data(length(first_half)+1:end)/ns;
          end
        end
        new_ads{j} = cur_ads;
      end
      perm_ads{k} = new_ads;
    end
  end
  
  
  
