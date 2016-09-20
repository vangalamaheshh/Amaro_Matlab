function scores = score_ziggs(Q,P,noise_thresh,broad_cutoff,length_hists)
  
  scores = zeros(size(Q,1),1);
    
  for i=1:size(Q,1)
    if Q(i,4) > 0 %% amplification
      sign = 1;
      if Q(i,4) > noise_thresh
        j = 1; %% significant amp
      else
        j = 2; %% insignificant amp
      end
    else %% deletion
      sign = -1;
      if Q(i,4) > -1*noise_thresh
        j = 3; %% insignificant del
      else
        j = 4; %% significant del
      end
    end
    scores(i) = P(j,1);
    
    % Now test if broad or focal!
    if Q(i,8) >= broad_cutoff
      scores(i) = scores(i) + P(j,2) + P(j,5) - P(j,4)*Q(i,4)*sign;
      lhist = length_hists{2};
    else
      scores(i) = scores(i) + P(j,3) + P(j,7) - P(j,6)*Q(i,4)*sign;
      lhist = length_hists{1};
    end
    
    ind = min(find(lhist{1} > Q(i,8)));
    
    if isempty(ind)
      ind = length(lhist{1});
    end
    
    scores(i) = scores(i) + log(lhist{2}(ind));
    
  end
  
      
    
