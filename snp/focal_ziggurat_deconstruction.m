function [QA QD QAOD QDOA] = focal_ziggurat_deconstruction(Amp,sample_broad_levels,with_cent,centromere_snps)
 
  cur_level = 0;
  unique_chrn = unique(Amp(:,1));
  for i=1:length(unique_chrn)
    ch = unique_chrn(i);
    ch_idx = find(with_cent == ch);
    if ~isempty(ch_idx) %% This means this chromosome has a centromere!
      idx = find(Amp(:,1)==ch);
      cent_idx = idx(find(Amp(idx,3)==centromere_snps(ch_idx)));
      %% p-arm
      cur_level = cur_level+1;
      cur_broad = sample_broad_levels(cur_level);
      AmpB = Amp(idx:cent_idx,:);
      [Qamp{cur_level} Qdel{cur_level} Qaod{cur_level} Qdoa{cur_level}] = ...
          zigg_deconstruction_loop(AmpB,cur_broad);
      
      %% q-arm
      cur_level = cur_level+1;
      cur_broad = sample_broad_levels(cur_level);
      AmpB = Amp((cent_idx+1):idx(end),:);
      [Qamp{cur_level} Qdel{cur_level} Qaod{cur_level} Qdoa{cur_level}] = ...
          zigg_deconstruction_loop(AmpB,cur_broad);
    else % This chromosome does not have a centromere!
      idx = find(Amp(:,1)==ch);
      cur_level = cur_level+1;
      cur_broad = sample_broad_levels(cur_level);
      AmpB = Amp(idx,:);
      [Qamp{cur_level} Qdel{cur_level} Qaod{cur_level} Qdoa{cur_level}] = ...
          zigg_deconstruction_loop(AmpB,cur_broad);
    end
    
  end

  QA = cat(1,Qamp{:});
  QD = cat(1,Qdel{:});
  QAOD = cat(1,Qaod{:});
  QDOA = cat(1,Qdoa{:});
