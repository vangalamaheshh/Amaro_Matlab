function Ampzig = focal_zigg_deconstruction(Amp)
  
  sample = Amp(1,5);
  Ampzig = [];
  Amptemp = [Amp(:,1:4) Amp(:,6)];
  Amptemp = merge_adj_segs(Amptemp);
  while length(Amptemp(:,4)) > 1 && max(Amptemp(:,4)) > 0
    [m,mi] = max(Amptemp(:,4));
    if mi == 1 %% max is leftmost segment on chromosome
      diff = Amptemp(mi,4)-Amptemp(mi+1,4);
      if diff>0
        Ampzig(size(Ampzig,1)+1,:) = [Amptemp(mi,1:3) diff sample Amptemp(mi+1,4) ...
                   Amptemp(mi,4:5)];
      end
      Amptemp(mi,:) = [Amptemp(mi,1:2) Amptemp(mi+1,3:4) Amptemp(mi,5)+Amptemp(mi+1,5)];
      Amptemp=Amptemp(setdiff(1:size(Amptemp,1),mi+1),:);
    elseif mi == size(Amptemp,1) %% max is rightmost segment on
                                 %chromosome
      diff = Amptemp(mi,4)-Amptemp(mi-1,4);
      if diff>0
      Ampzig(size(Ampzig,1)+1,:) = [Amptemp(mi,1:3) diff sample Amptemp(mi-1,4) ...
                   Amptemp(mi,4:5)]; 
      end
      Amptemp(mi-1,:) = [Amptemp(mi-1,1:2) Amptemp(mi,3) Amptemp(mi-1,4) Amptemp(mi,5)+Amptemp(mi-1,5)];
      Amptemp = Amptemp(1:mi-1,:);
    else %% max is in between other segments on chromosome
      if Amptemp(mi-1,4) >= Amptemp(mi+1,4)
        diff = Amptemp(mi,4)-Amptemp(mi-1,4);
        adj_idx = mi-1;
      else
        diff = Amptemp(mi,4)-Amptemp(mi+1,4);
        adj_idx = mi+1;
      end
      if diff>0
        Ampzig(size(Ampzig,1)+1,:)= [Amptemp(mi,1:3) diff sample Amptemp(adj_idx,4) ...
                          Amptemp(mi,4:5)];
      end
      if adj_idx < mi
        Amptemp(adj_idx,:) = [Amptemp(adj_idx,1:2) Amptemp(mi,3) ...
                            Amptemp(adj_idx,4) Amptemp(mi,5)+Amptemp(adj_idx,5)];
        Amptemp = Amptemp(setdiff(1:size(Amptemp,1),mi),:);
      else
        Amptemp(mi,:) = [Amptemp(mi,1:2) Amptemp(adj_idx,3:4) Amptemp(mi,5)+Amptemp(adj_idx,5)];
        Amptemp = Amptemp(setdiff(1:size(Amptemp,1),adj_idx),:);
      end
    end
  end
  %Amptemp = merge_adj_segs(Amptemp);
  if Amptemp(1,4)~=0
    Ampzig(size(Ampzig,1)+1,:) = [Amptemp(1,1:4) sample 0 Amptemp(1,4:5)];
  end
  
   
