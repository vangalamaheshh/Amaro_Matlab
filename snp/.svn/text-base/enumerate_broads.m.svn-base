function [broad_levels breakpoints] = enumerate_broads(D,cyto)
  
%% P = 4 x 7 matrix where 
%%  row 1 = scores for amps
%%  row 2 = scores for insig amp
%%  row 3 = scores for insig del 
%%  row 4 = scores for dels
%%  For each row:
%%    column 1 = event score
%%    column 2/3 = broad/focal score
%%    column 4/5 = alpha/beta_broad
%%    column 6/7 = alpha/beta_focal
%%    
  
  P = zeros(4,7);
  P(:,1)=[log(.15); log(.7); log(.7); log(.15)]; %%
                                                %amp/insig_amp/insig_del/del score
  P(:,2) = [log(.13); log(.5); log(.5); log(.13)]; %% broad_score
  P(:,3) = [log(1-.13); log(.5); log(.5);log(1-.13)]; %% focal scores
  P(:,4) = [4.3 45.4 43.2 5.45] %% broad alpha's
  P(:,5) = [log(6.6) log(22.4) log(22.4) log(9.4)]; %% log_beta broads
  P(:,6) = [4.66 18.2 17.6 4.9]; %% focal alphas
  P(:,7) = [log(7.43) log(10.7) log(10.7) log(8)]; %% log_beta_focals
   
  noise_thresh = .1;
  broad_cutoff = 0.99;
    
  [QA QD] = ziggurat(D,cyto);
  [fract_chr_armA chrarms] = normalize_by_arm_length(D,QA,cyto,1,2);
  fract_chr_armD = normalize_by_arm_length(D,QD,cyto,1,2,chrarms);
  
  [na,xout] = hist(fract_chr_armA,100);
  [nd,xout] = hist(fract_chr_armD,100);
  ntot = na+nd;
  
  focal_lengths = xout(1:49);
  broad_lengths = xout(50:end);
  
  focal_length_hist = ntot(1:49)/sum(ntot(1:49));
  broad_length_hist = ntot(50:end)/sum(ntot(50:end));
  
  length_hists{1} = {focal_lengths,focal_length_hist};
  length_hists{2} = {broad_lengths,broad_length_hist};
    
  broad_levels = zeros(2*max(D.chrn),size(D.dat,2));
  D.dat = 2.^(D.dat+1)-2;
  Qs = cell(1,size(D.dat,2));
  chr_bpts = zeros(max(D.chrn),size(D.dat,2));
  for ch=1:max(D.chrn) %% ch representes each chromosome number
    disp(['Finding broads on chromosome ' num2str(ch)])
    chr_snps = find(D.chrn == ch);
    %cent_snp = centromeres(ch);
    for j=1:size(D.dat,2) %% j represents each sample
      if mod(j,100)==0
        disp(j)
      end
      cur_chr = D.dat(chr_snps,j);
      bpts = find(diff(cur_chr)~=0)+min(chr_snps);
      bpts = union(bpts,max(chr_snps));
      B=zeros(length(bpts),6);
      B(:,1) = repmat(ch,length(bpts),1);
      B(:,5) = repmat(j,length(bpts),1);
      B(1,2) = min(chr_snps);
      B(2:end,2) = bpts(1:end-1)+1;
      B(:,3)=bpts;
      B(:,4)=D.dat(B(:,2),j);
      B(:,6) = normalize_by_arm_length(D,B,cyto,1,2,chrarms);
      
      %% FIXME: loop on each breakpoint
      
      if length(bpts) > 1
         bpt_scores = zeros(1,length(bpts));
         zigg_scores = zeros(1,length(bpts));
         pen_scores = zeros(1,length(bpts));
         p_levels = zeros(1,length(bpts));
         q_levels = zeros(1,length(bpts));
         max_Qs = cell(1,length(bpts));
         for i=1:length(bpts)
            if mod(i,10) == 0
              disp([j i])
            end
            bpt_snp = find(B(:,3)==bpts(i));
            if bpts(bpt_snp) == max(chr_snps)
              BP = B;
              p_fract = 2;
              BQ = [];
            else
              BP = B(1:bpt_snp,:);
              p_fract = sum(BP(:,6));
              BQ = B(bpt_snp+1:end,:);
            end
          
            if p_fract < broad_cutoff %% only q-arm is broad
                                      %disp('A');
              %[q_levels(i) max_Q q_score num_levels] = find_max_broad_level(D,BQ,P, ...
              %                                                  noise_thresh, ...
              %                                                  broad_cutoff, ...
              %                                                  length_hists, ...
              %                                                  2-p_fract);
              [q_levels(i) max_Q q_score num_levels] = find_max_broad_level_by_table(D,BQ,hd,xamp,ylen,2-p_fract);
              
              max_Q(:,10) = repmat(q_levels(i),size(max_Q,1),1);
              
              [QAp QDp] = ziggurat_on_extremes(D,BP,0,0,hd,xamp,ylen);
              
              if ~isempty(QAp)
                QAp = QAp(find(QAp(:,4)>=0),:);
              end
              
              if ~isempty(QDp)
                QDp = QDp(find(QDp(:,4)<0),:);
              end
              
              Qp = cat(1,QAp,QDp);
              
              Qp(:,10) = repmat(0,size(Qp,1),1);
              
              if ~isempty(Qp)
                p_score = sum(Qp(:,9));
              else
                p_score = 0;
              end
              
              zigg_scores(i) = p_score+q_score;
              pen_scores(i) = log(num_levels)+log(length(bpts)-1);
              %pen_scores(i) = log(num_levels);
              bpt_scores(i)=zigg_scores(i)-pen_scores(i);
              p_levels(i) = 0;
              max_Qs{i} = cat(1,Qp,max_Q);
              
            elseif p_fract >= broad_cutoff && p_fract <= (2-broad_cutoff)
              %disp('B');
              %% we can have a broad on either p or q-arms
              %[p_levels(i) max_Qp p_score num_levels_p] = find_max_broad_level(D,BP,P, ...
              %                                                  noise_thresh, ...
              %                                                  broad_cutoff, ...
              %                                                  length_hists, ...
              %                                                  p_fract);
              [p_levels(i) max_Qp p_score num_levels_p] = find_max_broad_level_by_table(D,BP,hd,xamp,ylen,p_fract);

              max_Qp(:,10) = repmat(p_levels(i),size(max_Qp,1),1);
              %[q_levels(i) max_Qq q_score num_levels_q] = find_max_broad_level(D,BQ,P, ...
              %                                                  noise_thresh, ...
              %                                                  broad_cutoff, ...
              %                                                  length_hists, ...
              %                                                  2- ...
              %                                                  p_fract);
              [q_levels(i) max_Qq q_score num_levels_q] = find_max_broad_level_by_table(D,BQ,hd,xamp,ylen,2-p_fract); 
 
              max_Qq(:,10) = repmat(q_levels(i),size(max_Qq,1),1);
 
              zigg_scores(i) = p_score+q_score;
              pen_scores(i) = log(num_levels_p)+log(num_levels_q);
              bpt_scores(i) = zigg_scores(i)-pen_scores(i);
              max_Qs{i} = cat(1,max_Qp,max_Qq);
              
            elseif p_fract > broad_cutoff && lt(p_fract,2) %% only p-arm is
                                                           %broad
                                                           %disp('C');
              %[p_levels(i) max_Q p_score num_levels] = find_max_broad_level(D,BP,P, ...
              %                                                  noise_thresh, ...
              %                                                  broad_cutoff, ...
              %                                                  length_hists, ...
              %                                                  p_fract);
              
              [p_levels(i) max_Q p_score num_levels] = find_max_broad_level_by_table(D,BP,hd,xamp,ylen,p_fract);


              max_Q(:,10) = repmat(p_levels(i),size(max_Q,1),1);
              [QAq QDq] = ziggurat_on_extremes(D,BQ,0,0,hd,xamp,ylen);
              
              if ~isempty(QAq)
                QAq = QAq(find(QAq(:,4)>0),:);
              end
              if ~isempty(QDq)
                QDq = QDq(find(QDq(:,4)<0),:);
              end
              
              Qq = cat(1,QAq,QDq);
              
              Qq(:,10) = repmat(0,size(Qq,1),1);
              
              if ~isempty(Qq)
                q_score = sum(Qq(:,9));
              else
                q_score = 0;
              end
              zigg_scores(i) = p_score+q_score;
              pen_scores(i) = log(num_levels)+log(length(bpts)-1);
              %pen_scores(i) = log(num_levels);
              bpt_scores(i)=zigg_scores(i)-pen_scores(i);
              q_levels(i) = 0;
              max_Qs{i} = cat(1,max_Q,Qq);
            else %% p-fract == 2, and bpt is at telomere
                 %disp('D');
              %[p_levels(i) max_Qs{i} max_score num_levels] = find_max_broad_level(D,BP,P, ...
              %                                                  noise_thresh, ...
              %                                                  broad_cutoff, ...
              %                                                  length_hists, ...
              %                                                  p_fract);
              [p_levels(i) max_Qs{i} max_score num_levels] = find_max_broad_level_by_table(D,BP,hd,xamp,ylen,p_fract);      
              max_Qs{i}(:,10) = repmat(p_levels(i),size(max_Qs,1),1);
              zigg_scores(i) = max_score;
              pen_scores(i)=log(num_levels);
              bpt_scores(i) = zigg_scores(i)-pen_scores(i);
              q_levels(i) = p_levels(i);
            end
            [mx mk1] = max(bpt_scores);
            [mx mk] = max(zigg_scores);
            Qs{j} = cat(1,Qs{j},max_Qs{mk});
            chr_bpts(ch,j) = bpts(mk);
            broad_levels(2*ch-1,j) = p_levels(mk);
            broad_levels(2*ch,j) = q_levels(mk);
         end
      else
        disp('Only 1 segment on chromosome!')
        Qs{j} = [B(1,1:5) 0 B(1,4) 2 0 B(1,4)];
        chr_bpts(ch,j) = max(chr_snps);
        broad_levels(2*ch-1,j) = B(1,4);
        broad_levels(2*ch,j) = B(1,4);
        mk=1;
        mk1 = 1;
      end
      plot_ziggs(cur_chr,Qs{j},broad_levels(2*ch-1,j),broad_levels(2*ch, ...
                                                        j),chr_snps(1),bpts(mk),chr_snps(end));
      if mk ~= mk1
        disp('Penalty term in play!')
      end
      keyboard
    end
  end
  
