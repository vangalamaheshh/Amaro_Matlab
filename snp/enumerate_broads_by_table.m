function [D QA QD QAOD QDOA hd2] = enumerate_broads_by_table(D,cyto,do_plot,max_segs_per_sample)
  
  if ~exist('do_plot','var') || isempty(do_plot)
    do_plot = 0;
  end
  
  if ~exist('max_segs_per_sample','var') || isempty('max_segs_per_sample')
    max_segs_per_sample = 500;
  end
  
%% Ziggurat deconstruction (w respect to 0 level)
  
  [QA QD num_bpts] = ziggurat(D,cyto);
  QD(:,4) = -1*QD(:,4);
  Q = cat(1,QA,QD);
  
  %% Normalize by arm length
  [fract_chr_armA chrarms] = normalize_by_arm_length(D,QA,cyto,1,2);
  fract_chr_armD = normalize_by_arm_length(D,QD,cyto,1,2,chrarms);
  fract = cat(1,fract_chr_armA,fract_chr_armD);
  
  %% Generate 2-D length/amplitude histogram
  xamp = -2:.08:2;
  ylen = 0:0.04:2;
  
  %good = find(fract > .05);
  
  %hd = hist2d(Q(good,4),fract(good),xamp,ylen);     
  hd = hist2d(Q(:,4),fract,xamp,ylen);
  hd1 = hd+1;
  hd2 = hd1/sum(sum(hd1));     
  
  log_hd = log(hd2);
  
  %% Plot 2-D length/amplitude histogram
  
  bar3(log(hd2)+15)
  set(gca,'XTickLabel',{.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2})
  set(gca,'YTickLabel',{-2,-1.2,-.4,.4,1.2,2})
  xlabel('Length (fract of chr arm)')
  ylabel('Amplitude')
  zlabel('log(freq)');
  title('Distribution of segments as function of length and amplitude')
  
  %% Initialize variables
  
  broad_levels = zeros(2*max(D.chrn),size(D.dat,2)); %% holds broad levels for each arm 
  D.dat = 2.^(D.dat+1)-2; %% convert data to copy number space
  Qs = cell(max(D.chrn),size(D.dat,2)); %% holds max ziggurat for each sample and each chromosome 
  chr_bpts = zeros(max(D.chrn),size(D.dat,2)); %% snp location of
                                               %breakpoint on each
                                               %chromosome in each sample
  
  
  %% LEFT OFF HERE!
      
  %[bpt_ids sample_ids] = find(diff(D.dat) ~=0);
  
  %for i=1:size(D.dat,2)
  %  num_bpts(i) = length(find(sample_ids==i));
  %end
  
  samples_to_rm = find(num_bpts> max_segs_per_sample);
  
  D = reorder_D_cols(D,setdiff(1:size(D.dat,2),samples_to_rm));
  
  D.zigg_id = 1:size(D.dat,2);
  for ch=1:max(D.chrn) %% ch representes each chromosome number
    disp(['Finding broads on chromosome ' num2str(ch)])
    chr_snps = find(D.chrn == ch);
    for j=1:size(D.dat,2) %% j represents each sample
      if mod(j,100)==0
        disp(j)
      end
      %% get data from this chromosome, and initialize B for iterative ziggurat deconstruction
      cur_chr = D.dat(chr_snps,j);
      bpts = find(diff(cur_chr)~=0)+min(chr_snps); 
      %good = find(~isnan(D.dat(bpts,j)));
      %bpts = bpts(good);
      bpts = union(bpts,max(chr_snps));
      B=zeros(length(bpts),6);
      B(:,1) = repmat(ch,length(bpts),1);
      B(:,5) = repmat(j,length(bpts),1);
      B(1,2) = min(chr_snps);
      B(2:end,2) = bpts(1:end-1)+1;
      B(:,3)=bpts;
      B(:,4)=D.dat(B(:,2),j);
      B(:,6) = normalize_by_arm_length(D,B,cyto,1,2,chrarms);
                  
      %% loop over each breakpoint and find most likely broad level on
      %each arm for each breakpoint
      %% Then, find most likely breakpoint overall
      
      if length(bpts) > 1
        bpt_scores = zeros(1,length(bpts));
        zigg_scores = zeros(1,length(bpts));
        pen_scores = zeros(1,length(bpts));
        p_levels = zeros(1,length(bpts));
        q_levels = zeros(1,length(bpts));
        p_fract = zeros(1,length(bpts));
        q_fract = zeros(1,length(bpts));
        max_Qs = cell(1,length(bpts));
        for i=1:length(bpts)
          if mod(i,10) == 0
            disp([j i])
          end
          bpt_snp = find(B(:,3)==bpts(i));
          if bpts(bpt_snp) == max(chr_snps)
            BP = B;
            p_fract(i) = sum(B(:,6));
            BQ = [];
            q_fract(i) = 0;
          else
            BP = B(1:bpt_snp,:);
            p_fract(i) = sum(BP(:,6));
            BQ = B(bpt_snp+1:end,:);
            q_fract(i) = sum(BQ(:,6));
          end
          
          [p_levels(i) max_Qp p_score num_levels_p] = ...
              find_max_broad_level_by_table(D,BP,log_hd,xamp,ylen,p_fract(i));
          
          max_Qp(:,10) = repmat(p_levels(i),size(max_Qp,1),1);
          
          if ~isempty(BQ)
            [q_levels(i) max_Qq q_score num_levels_q] = find_max_broad_level_by_table(D,BQ,log_hd,xamp,ylen,q_fract(i));
            max_Qq(:,10) = repmat(q_levels(i),size(max_Qq,1),1);
            len_bpts = length(bpts);
          else
            q_levels(i) = p_levels(i);
            max_Qq = [];
            q_score = 0;
            num_levels_q = 1;
            len_bpts = 1;
          end
          
          zigg_scores(i) = p_score+q_score;
          pen_scores(i) = log(num_levels_p+1) + log(num_levels_q+1)+ ...
              log(len_bpts);
          bpt_scores(i)=zigg_scores(i)-pen_scores(i);
          max_Qs{i} = cat(1,max_Qp,max_Qq);
          
        end
        [mx mk] = max(bpt_scores);
        [mx1 mk1] = max(zigg_scores);
        Qs{ch,j} = max_Qs{mk};
        chr_bpts(ch,j) = bpts(mk);
        broad_levels(2*ch-1,j) = p_levels(mk);
        broad_levels(2*ch,j) = q_levels(mk);
      else
        disp('Only 1 segment on chromosome!')
        Qs{ch,j} = [B(1,1:5) 0 B(1,4) sum(B(:,6)) 0 B(1,4)];
        chr_bpts(ch,j) = max(chr_snps);
        broad_levels(2*ch-1,j) = B(1,4);
        broad_levels(2*ch,j) = B(1,4);
        mk=1;
        mk1 = 1;
      end
        if do_plot
          plot_ziggs(cur_chr,Qs{ch,j},broad_levels(2*ch-1,j),broad_levels(2*ch, ...
                                                            j),chr_snps(1),bpts(mk),chr_snps(end));
          keyboard
        end
    end
    
  end
  
for ch=1:max(D.chrn)
  Qs_temp{ch} = cat(1,Qs{ch,:});
end
  
Q = cat(1,Qs_temp{:});
  
QA = Q(find(Q(:,4) > 0 & Q(:,6) >= 0),:);
QD = Q(find(Q(:,4) < 0 & Q(:,6) <= 0),:);

QAOD = Q(find(Q(:,4) > 0 & Q(:,7)<=0),:);
QDOA = Q(find(Q(:,4) < 0 & Q(:,7)>=0),:);

spans_zeroA = find(Q(:,4) >0 & Q(:,6).*Q(:,7) < 0);
spans_zeroD = find(Q(:,4) < 0 & Q(:,6).*Q(:,7)<0);

new_a = cell(1,length(spans_zeroA));
new_aod = cell(1,length(spans_zeroA));
for i=1:length(spans_zeroA)
  cur_seg = Q(spans_zeroA(i),:);
  aod_seg = [cur_seg(1:6) 0 cur_seg(8:10)];
  amp_seg = [cur_seg(1:5) 0 cur_seg(7:10)];
  aod_seg(4) = aod_seg(7)-aod_seg(6);
  amp_seg(4) = amp_seg(7)-amp_seg(6);
  new_a{i} = amp_seg;
  new_aod{i} = aod_seg;
end

new_d = cell(1,length(spans_zeroD));
new_doa = cell(1,length(spans_zeroD));
for i=1:length(spans_zeroD)
  cur_seg = Q(spans_zeroD(i),:);
  doa_seg = [cur_seg(1:6) 0 cur_seg(8:10)];
  del_seg = [cur_seg(1:5) 0 cur_seg(7:10)];
  doa_seg(4) = doa_seg(7)-doa_seg(6);
  del_seg(4) = del_seg(7)-del_seg(6);
  new_d{i} = del_seg;
  new_doa{i} = doa_seg;
end

QA = cat(1,QA,new_a{:});
QD = cat(1,QD,new_d{:});

QAOD = cat(1,QAOD,new_aod{:});
QDOA = cat(1,QDOA,new_doa{:});
  



