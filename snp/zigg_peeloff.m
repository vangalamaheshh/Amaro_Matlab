function Ampzig = zigg_peeloff(Amp)

  
%% The zigg_peeloff function takes a matrix Amp which is assumed to contain
%% the (amp or del) segments for a single sample and performs iterative peel-off to generate a set of
%% segments reflecting the underlying genomic fragments giving rise to that 
%% pattern of segmental alteration
  
%% Amp(:,1) = chrn of segment
%% Amp(:,2) = starting snp of segment
%% Amp(:,3) = ending snp of segment
%% Amp(:,4) = log2(val) of segment  

  Ampzig = [];
  for ch=1:22
    idx = find(Amp(:,1)==ch);
    % check if empty
    Amptemp=Amp(idx,1:4);
    Amptemp = merge_adj_segs(Amptemp);
    while length(Amptemp(:,4)) > 1 && max(Amptemp(:,4)) > 0
      [m,mi] = max(Amptemp(:,4));
      if mi == 1 %% max is leftmost segment on chromosome
        diff = Amptemp(mi,4)-Amptemp(mi+1,4);
        Ampzig(size(Ampzig,1)+1,:) = [Amptemp(mi,1:3) diff];  
        Amptemp(mi,:) = [Amptemp(mi,1:2) Amptemp(mi+1,3:4)];
        Amptemp=Amptemp(setdiff(1:size(Amptemp,1),mi+1),:);
      elseif mi == size(Amptemp,1) %% max is rightmost segment on
                                    %chromosome
        diff = Amptemp(mi,4)-Amptemp(mi-1,4);
        Ampzig(size(Ampzig,1)+1,:) = [Amptemp(mi,1:3) diff];  
        Amptemp(mi-1,:) = [Amptemp(mi-1,1:2) Amptemp(mi,3) Amptemp(mi-1,4)];
        Amptemp = Amptemp(1:mi-1,:);
      else %% max is in between other segments on chromosome
        if Amptemp(mi-1,4) >= Amptemp(mi+1,4)
          diff = Amptemp(mi,4)-Amptemp(mi-1,4);
          adj_idx = mi-1;
        else
          diff = Amptemp(mi,4)-Amptemp(mi+1,4);
          adj_idx = mi+1;
        end
        Ampzig(size(Ampzig,1)+1,:)= [Amptemp(mi,1:3) diff];
        if adj_idx < mi
          Amptemp(adj_idx,:) = [Amptemp(adj_idx,1:2) Amptemp(mi,3) Amptemp(adj_idx,4)];
          Amptemp = Amptemp(setdiff(1:size(Amptemp,1),mi),:);
        else
          Amptemp(mi,:) = [Amptemp(mi,1:2) Amptemp(adj_idx,3:4)];
          Amptemp = Amptemp(setdiff(1:size(Amptemp,1),adj_idx),:);
        end
      end
    end
    Amptemp = merge_adj_segs(Amptemp);
    Ampzig(size(Ampzig,1)+1,:) = Amptemp(1,:);
  end
   
  sample_id = repmat(Amp(1,5),size(Ampzig,1),1);
  Ampzig = cat(2,Ampzig,sample_id);
  
