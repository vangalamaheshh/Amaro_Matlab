function [Qa Qd Qaod Qdoa] = zigg_deconstruction_loop(AmpB,cur_broad)
  
  Qa = [];
  Qd = [];
  Qaod = [];
  Qdoa = [];
  
  if cur_broad == 0
    DelB = AmpB;
    AmpB(:,4) = AmpB(:,4).*(AmpB(:,4)>0);
    DelB(:,4) = -1*DelB(:,4).*(DelB(:,4)<0);
    Qa = focal_zigg_deconstruction(AmpB);
    Qd = focal_zigg_deconstruction(DelB);
  elseif cur_broad > 0
    DelB = AmpB;
    AmpB(:,4) = AmpB(:,4).*(AmpB(:,4)>0); %% zeros everything below zero
    DelB(:,4) = -1*DelB(:,4).*(DelB(:,4)<0); %% below zero is a deletion
    Qd = focal_zigg_deconstruction(DelB);
    
    %% Now subtract median level and find residual amplifications and
    %DOAs
    AmpB(:,4) = AmpB(:,4)-repmat(cur_broad,size(AmpB,1),1);
    DoaB = AmpB;
    AmpB(:,4) = AmpB(:,4).*(AmpB(:,4)>0); %% zeros everything below zero
    DoaB(:,4) = -1*DoaB(:,4).*(DoaB(:,4)<0);
    Qa = focal_zigg_deconstruction(AmpB);
    Qdoa = focal_zigg_deconstruction(DoaB);
    
  else
    DelB = AmpB;
    AmpB(:,4) = AmpB(:,4).*(AmpB(:,4)>0); %% zeros everything below zero
    DelB(:,4) = -1*DelB(:,4).*(DelB(:,4)<0); %% below zero is a deletion
    Qa = focal_zigg_deconstruction(AmpB);
    
    %% Now subtract median level and find residual deletions and
    % AODs
    DelB(:,4) = DelB(:,4)-repmat(-1*cur_broad,size(DelB,1),1);
    AodB = DelB;
    DelB(:,4) = DelB(:,4).*(DelB(:,4)>0); %% zeros everything below zero
    AodB(:,4) = -1*AodB(:,4).*(AodB(:,4)<0);
    Qd = focal_zigg_deconstruction(DelB);
    Qaod = focal_zigg_deconstruction(AodB);
  end
  
  
