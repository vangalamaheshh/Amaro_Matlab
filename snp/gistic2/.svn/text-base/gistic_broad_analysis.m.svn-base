function [arm_medians qA qD zA zD fA fD names num_genes params] = ...
                gistic_broad_analysis(base_dir,D,ref_gene_file,params)
%GISTIC_BROAD_ANALYSIS - analyze arm-level events

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% Set default parameters
  params = impose_default_value(params,'t_amp',0.1);
  params = impose_default_value(params,'t_del',0.1);
  params = impose_default_value(params,'ext','');
  params = impose_default_value(params,'fname','');
  params = impose_default_value(params,'use_two_sided',0);
  params = impose_default_value(params,'broad_len_cutoff',0.98);
  params = impose_default_value(params,'genepattern',0);
  thresh = [params.t_amp params.t_del params.t_amp params.t_del];

  use_segarray = isa(D.dat,'SegArray');

  % Initalize variables
  n = size(D.dat,2);
  load(ref_gene_file);
  
  fname = params.fname;
  ext = params.ext;
  
  broads = reconstruct_genomes(D.Qs,struct('broad_or_focal','broad',...
                                           'broad_len_cutoff',params.broad_len_cutoff,...
                                           't_amp',params.t_amp,...
                                           't_del',params.t_del,...
                                           'column_to_add',12,...
                                           'use_segarray',use_segarray,...
                                           'rows',length(D.pos)));

  %% Amplifications
  
  Z.dat = broads.amp;
  Z.chrn = D.chrn; Z.pos = D.pos;
  
  verbose('Calculating median of arm values...',20);
  % find medians and retain those over the significance threshold
  arm_mediansA = find_med_chrarm(Z,cyto);
  arm_mediansA(arm_mediansA < thresh(1)) = 0;
  
  % find arms with enough markers
  w_data = find(~isnan(arm_mediansA(:,1)));
  w_data = setdiff(w_data,41); %% get rid of 21p, which has only 8 genes !!!HUMCHR
                               %and few snps
  
  %% Deletions
  
  Z.dat = broads.del;
  arm_mediansD = find_med_chrarm(Z,cyto);
  arm_mediansD(arm_mediansD < thresh(1)) = 0;
  
  % Calculate arm_scores, frequency
  arm_medians = arm_mediansA - arm_mediansD;
  arm_medians = arm_medians(w_data,:);
  
  arm_mediansA(arm_mediansD > 0) = NaN;
  arm_mediansD(arm_mediansA > 0) = NaN;
  
  % Arm frequencies
  fA = sum(arm_mediansA(w_data,:) > 0,2)./(n-sum(arm_mediansD(w_data,:) > 0,2));
  fD = sum(arm_mediansD(w_data,:) > 0,2)./(n-sum(arm_mediansA(w_data,:) > ...
                                                 0,2));
    
  %% Calculate number of genes on each chromosome arm

  armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9XY]+[pq]+','match'),'UniformOutput',false));
  armnames = armnames(~cellfun(@isempty,armnames));
  
  chrarms = struct('name',{},'start',{},'end',{},'chrn',{},'length',{});
  
  for i=1:length(armnames)
    idx = strmatch(armnames(i),{cyto.name});
    band.name = char(armnames(i));
    band.start = cyto(idx(1)).start+1;
    band.end = cyto(idx(end)).end;
    band.chrn=cyto(idx(end)).chrn;
    band.length = band.end-band.start+1;
    chrarms(i) = band;
  end
  
  chr=cat(1,chrarms.chrn);
  [~,sposi]=sort(chr);
  chrarms = chrarms(sposi);
  
  arm_data = chrarms(w_data);
  num_genes = zeros(1,length(arm_data));
  names = {arm_data.name};
  
  for i=1:length(arm_data)
    disp(i)
    num_genes(i) = length(find([rg.chrn] == arm_data(i).chrn & [rg.start] >= ...
                               arm_data(i).start & [rg.end] <= ...
                               arm_data(i).end));
  end
  
  %% scatter plot of number-of-genes vs frequency-of-amplification
  h = figure();
  tufte_plots(num_genes,(fA+fD)/2,names,h);
  title('Correlation of frequency of arm-level events vs. number of genes on each chromosome arm');
  xlabel('Number of genes on chr. arm');
  ylabel('Frequency of arm-level events');
  
  freqy_name = [base_dir fname 'freqarms_vs_ngenes' ext];
  if ~params.genepattern
      saveas(gcf,[freqy_name '.fig'],'fig');
  end
  saveas(gcf,[freqy_name '.pdf'],'pdf');
  
  %% calculate z-scores and p/q-values
  
  b = robustfit(num_genes',.5*(fA+fD));
  % calculate expected frequencies and make sure we're powered
  powered = true;
  f_hat = b(2)*num_genes'+b(1);
  if any(f_hat <= 0)
    warning('Too few events to power GISTIC broad significance analysis');
    f_hat = max(0,f_hat);
    powered = false;
  end
  
  % get count of all arms with medians
  nsA = sum(~isnan(arm_mediansA(w_data,:)),2);
  nsD = sum(~isnan(arm_mediansD(w_data,:)),2);
  
  xA = nsA.*fA;
  xD = nsD.*fD;
  
  % predicted mean number of alterations based on fitted line
  uA = nsA.*f_hat;
  uD = nsD.*f_hat;
  
  % predicted standard deviations based on binomial distribution
  sA = sqrt(nsA.*f_hat.*(1-f_hat));
  sD = sqrt(nsD.*f_hat.*(1-f_hat));
  
  % z-score of observed values
  zA = (xA-uA)./sA;
  zD = (xD-uD)./sD;
  
  if params.use_two_sided
    % convert to two-sided test
    pAu = 1-normcdf(zA,0,1);
    pDu = 1-normcdf(zD,0,1);
    
    pAd = normcdf(zA,0,1);
    pDd = normcdf(zD,0,1);
    
    [pA pAidx] = min(2*[pAu pAd],[],2);
    pAidx = 1-2*(pAidx-1);
    
    [pD pDidx] = min(2*[pDu pDd],[],2);
    pDidx = 1-2*(pDidx-1);
  else
    % one-sided test
    pA = 1-normcdf(zA,0,1);
    pD = 1-normcdf(zD,0,1);
  end
  % calculate BH false discovery rates
  qA = calc_fdr_value(pA);
  qD = calc_fdr_value(pD);
  
  %% for two-sided test, create amps versus dels scatter plot
  
  if params.use_two_sided && powered %% generate 2d quadrant figure
    % Find thresholds
    amp_pos_thresh = min(zA(intersect(find(pAidx == 1),find(qA <= .25))));
    amp_neg_thresh = max(zA(intersect(find(pAidx == -1),find(qA <= .25))));
    del_pos_thresh = min(zD(intersect(find(pDidx == 1),find(qD <= .25))));
    del_neg_thresh = max(zD(intersect(find(pDidx == -1),find(qD <= .25))));
    
    figure()
    w=2;
    f=10;
    dz=.25;
    scatter(zA,zD,'Marker','none')
    
    ymax = 10*(1+floor(max(abs(zD))/10));
    ymin = -1*ymax;
    xmax = 10*(1+floor(max(abs(zA))/10));
    xmin = -1*xmax;
    
    line([amp_pos_thresh amp_pos_thresh],[ymin 1*ymax],'Color','Green','LineStyle',':','LineWidth',w)
    line([amp_neg_thresh amp_neg_thresh],[ymin ymax],'Color','Green','LineStyle',':','lineWidth',w)
    line([xmin xmax],[del_pos_thresh del_pos_thresh],'Color','Green','LineStyle',':','LineWidth',w)
    line([xmin xmax],[del_neg_thresh del_neg_thresh],'Color','Green','LineStyle',':','LineWidth',w)
    
    for j=1:length(names)
      if zD(j) >= del_pos_thresh && zA(j) >= amp_pos_thresh
        text(zA(j),zD(j)+dz,names(j),'Color','magenta','FontSize',f);
      elseif zA(j) >= amp_pos_thresh 
        text(zA(j),zD(j)+dz,names(j),'Color','red','FontSize',f);
      elseif zD(j) >= del_pos_thresh 
        text(zA(j),zD(j)+dz,names(j),'Color','blue','FontSize',f);
      else
        text(zA(j),zD(j)+dz,names(j),'Color','black','FontSize',f);
      end
    end
    
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    
    xlabel('Amplification Z score')
    ylabel('Deletion Z score')
    
    % save image files
    if ~params.genepattern
      saveas(gcf,[base_dir fname 'broad_gistic_plot' ext '.fig'],'fig')
    end
    saveas(gcf,[base_dir fname 'broad_gistic_plot' ext '.pdf'],'pdf')
  end
  
  %% save broad results table
  
  f = fopen([base_dir fname 'broad_significance_results' ext '.txt'],'w');
  % issue low-frequency low-power warning
  if ~powered
      fprintf(f,'# WARNING: under-powered analysis due to low broad event frequency\n');
  end
  fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Arm','# Genes',...
          'Amp frequency','Amp z-score','Amp q-value',...
          'Del Frequency','Del z-score','Del q-value');
  
  for j=1:length(names)
    fprintf(f,'%s\t%d\t%0.2f\t%1.3g\t%1.3g\t%0.2f\t%1.3g\t%1.3g\n',char(names(j)), ...
            num_genes(j),fA(j),zA(j),qA(j),fD(j),zD(j),qD(j));
  end
  fclose(f);
  
  % also save as matlab file
  if ~params.genepattern
    save([base_dir fname 'broad_results' ext '.mat'],'zA','zD','fA','fD','qA', ...
       'qD','num_genes','names');
  end
  

    
