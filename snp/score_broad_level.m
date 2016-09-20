function broad_score = score_broad_level(broad_level,P,noise_thresh,prob_zero)
  
  %zeroscore = log(prob_zero);
  %ampscore = log((1-prob_zero)*exp(P(1,1)));
  ampscore = P(1,1);
  %delscore = log((1-prob_zero)*exp(P(2,1)));
  delscore = P(2,1);
  insig_beta_score = log(P(1,10)/(1-exp(-P(1,10)*noise_thresh)));
  
  if broad_level >0 % broad is an amplification
    broad_score = ampscore;
      if broad_level > noise_thresh
        broad_score = broad_score+P(1,2);
        broad_score = broad_score+P(1,4)+log(P(1,8))+P(1,8)*(noise_thresh-broad_level);
      else
        broad_score = broad_score+P(1,3)+P(1,6)+insig_beta_score-P(1,10)*broad_level;
      end
  elseif broad_level < 0 %% broad is a deletion
    broad_score = delscore;
    if broad_level < -1*noise_thresh %% broad deletion is significant
      broad_score = broad_score+P(2,2);
      broad_score = broad_score+P(2,4)+log(P(2,8))+P(2,8)*(noise_thresh+broad_level);
    else %% broad deletion is insignificant
      broad_score = broad_score+P(2,3)+P(2,6)+insig_beta_score+P(2,10)*broad_level;
    end 
  else %% no event
    %broad_score = zeroscore+insig_beta_score;
    broad_score = insig_beta_score;
  end
    
