function [D QA QD QAOD QDOA] = perform_deconstruction2(D,cyto,unused,niters)
%PERFORM_DECONSTRUCTION inner, low-level ziggurate deconstruction function
%
%   [D QA QD QAOD QDOA] = perform_deconstruction(D,CYTO,unused,niters)
%
% PERFORM_DECONSTRUCTION takes a D structure and cytoband information in 
% CYTO as inputs and performs ziggurate deconstruction on the data 
% using the number of iterations specified in NITERS. It is assumed that
% the CN data in D.dat are in (absolute copy number - 2) units.
%
% The results are copy number events returned in four Q arrays: QA for pure
% amplification events, QD for pure deletion events, QAOD for
% amplifications over deletions and QDOA for deletions over amplifications.
%
% see also Qs.
%
% The third argument to PERFORM_DECONSTRUCTION is no longer used and should
% be left empty.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%% General Steps: 1) Starting with arm_levels = 0, perform deconstruction
%%                2) Using data, find most likely arm-levels and centromeres
%%                3) Repeat
  
  
  %% Step 1: Given broad level (init 0) and bpt (init telomere), do deconstruction 
  %% Step 2: Calculate segment length-amp distribution
  %% Step 3: Given distribution, find broad level and bpt that maximizes
  % likelihood function
  %% Step 4: Repeat
  
  %% Need a function which takes a list of broad levels and potential
  %breakpoints, loops over them, and returns all deconstructions with
  %those parameters (initially, run with no broad level and telomeric bpt).
  
  %% Then, need to score all deconstructions given current
  %bg_distribution, and pick broad levels and breakpoints that maximizes
  %the likelihood of the data
  
  %% Then regenerate background model, and repeat (if necessary)

  if ~exist('niters','var') || isempty(niters)
    niters = 1;
  end
  
  % Initialize variables
  chrnEnd = [];
  chrarms = [];
  Zamp = cell(1,size(D.dat,2));
  Zdel = Zamp;
    
  %% Initial deconstruction (w respect to 0 level)
    
  broad_levels = zeros(2*max(D.chrn),size(D.dat,2));
    
  for j=1:size(D.dat,2)
    modi(j,100)
    
    % create data matrix
    [B chrnEnd chrarms] = make_sample_B(D,j,cyto,chrnEnd,chrarms);
    
    % do deconstruction
    [Zamp{j} Zdel{j}] = deconstruct_sample(B,broad_levels(:,j),chrnEnd);
  end
  
  QA = cat(1,Zamp{:}); 
  QD = cat(1,Zdel{:});
  QD(:,4) = -1*QD(:,4);
  clear Zamp Zdel
  
  %% Generate length-amp histogram
  xamp = [unique(QD(:,4));unique(QA(:,4))]';
  [log_hd xamp ylen] = generate_2d_hists(QA,QD,xamp,[],.01,1);
  log_hd_st = log_hd;  
  %% Loop over number of iterations
  
  for l=1:3
    Qs = cell(max(D.chrn),size(D.dat,2)); %% holds max ziggurat for each sample and each chromosome 
    chr_bpts = zeros(max(D.chrn),size(D.dat,2)); %% snp location of
                                                 %breakpoint on each
                                                 %chromosome in each sample
    
    
    % loop over samples
    for j=1:size(D.dat,2)
      modi(j,25)
      [B chrnEnd chrarms] = make_sample_B(D,j,cyto,chrnEnd,chrarms);
      % loop over chromosomes in each sample
      for ch=1:max(B(:,1))
        % get segments for current chromosome
        Bt = B(B(:,1) == ch,:);
        if size(Bt,1) == 1
          % Only 1 segment on chromosome
          Qs{ch,j} = [Bt(1,1:5) 0 Bt(1,4) sum(Bt(:,6)) 0 Bt(1,4) 0 0 Bt(1,13:14)];
          chr_bpts(ch,j) = Bt(1,3);
          broad_levels(2*ch-1,j) = Bt(1,4);
          broad_levels(2*ch,j) = Bt(1,4);
        else
          % initialize temporary variables
          bpt_scores = zeros(1,size(Bt,1));
          zigg_scores = bpt_scores;
          pen_scores = bpt_scores;
          p_levels = bpt_scores;
          q_levels = bpt_scores;
          p_fract = bpt_scores;
          q_fract = bpt_scores;
          max_Qs = cell(1,size(Bt,1));
          event_score = zeros(1,size(Bt,1));
          %% CHANGE: Adding event probability.  Assume Poisson with expected number
          % of events lamda.
          num_segments = size(Bt,1);
          % Estimating lamda to be expected number of events on the chromosome
          % given the segmentation.  Coarse estimate here: Each event causes 2
          % breakpoints.  Generally, breakpoints won't align exactly.  Might add
          % a bit to lamda to consider possibility? Additionally, perhaps
          % not unlikely for endpoints to be part of two separate events,
          % so again, too simple...
          if Bt(1,4) == 0 && Bt(end,4) == 0  % Include neither end as breakpoint.
              lamda = ceil((num_segments - 1)/2); 
          elseif Bt(1,4) == 0 || Bt(end,4) == 0  % Include 1 end.
              lamda = ceil((num_segments)/2);
          else
              lamda = ceil((num_segments +1)/2); % Include both.
          end
          for i=1:size(Bt,1) %% loop over all breakpoints
            % Calculate p,q fractions
            p_fract(i) = sum(Bt(1:i,8));
            q_fract(i) = sum(Bt(i+1:end,8)); 
            
            % Find max p,q levels for this breakpoint
            [p_levels(i) max_Qp p_score num_levels_p] = find_max_broad_level_by_table2(Bt(1:i,:),log_hd,xamp,ylen,p_fract(i));
            [q_levels(i) max_Qq q_score num_levels_q] = find_max_broad_level_by_table2(Bt(i+1:end,:),log_hd,xamp,ylen,q_fract(i));
            
            max_Qp(:,10) = repmat(p_levels(i),size(max_Qp,1),1);
            max_Qq(:,10) = repmat(q_levels(i),size(max_Qq,1),1);
            
            zigg_scores(i) = p_score+q_score;
            
            if num_levels_q == 0
              len_bpts = 1;
            else
              len_bpts = size(Bt,1);
            end
            %% Change: Adjusting 'penalty':
            % Add BIC penalty score = (2*k-1) * ln n where k= #broad levels
%             if num_levels_p == 0 || num_levels_q == 0
%               pen_scores(i) = log(len_bpts); % k = 1, so pen = ln(n)
%             else
%               pen_scores(i) = 3*log(len_bpts); % k = 2, so pen = 3*ln(n)
%             end
%             bpt_scores(i)=zigg_scores(i)-pen_scores(i);
            
            max_Qs{i} = cat(1,max_Qp,max_Qq);
            %% Add event count:
            num_events = size(max_Qs{i},1);
            if (p_levels(i) ~= 0) && (p_levels(i) == q_levels(i))
                num_events = num_events - 1;  % Merge two broad events if at same level.  
            end
            %% To do: adjust above merge in Qs before combining max_Qp and max_Qq so
            % that histogram will be more accurate.
            event_score(i) = num_events*log(lamda) - log(factorial(num_events));
            %Note: log-probabiliy = log((lamda*t)^num_events*exp(-lamda*t)/factorial(num_events))
            % comes from Poisson density taking t = 1.
            % Dropping '-lamda' term as identical for each decomposition.
            bpt_scores(i) = zigg_scores(i) + event_score(i);
            
          end
          [mx mk] = max(bpt_scores);
          Qs{ch,j} = max_Qs{mk};
          chr_bpts(ch,j) = Bt(mk,3);
          broad_levels(2*ch-1,j) = p_levels(mk);
          broad_levels(2*ch,j) = q_levels(mk);
        end
        
      end
      
    end
    
    
    for ch=1:max(D.chrn)
      Qs_temp{ch} = cat(1,Qs{ch,:});
    end
    
    Q = cat(1,Qs_temp{:});
    
    [QA QD QAOD QDOA] = make_final_Qs(Q);
    
    [log_hd xamp ylen] = generate_2d_hists(QA,QD,xamp,[],.01,1);
    hds{l} = log_hd;
  end
  
  
 
