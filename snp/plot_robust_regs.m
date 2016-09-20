function plot_robust_regs(D,regs,ads,nsamples,snp_margin,show_plots,gene_gistic_rg)
  
  if ~exist('show_plots','var') || isempty(show_plots)
    show_plots = 0;
  end
  
  if ~exist('nsamples','var') || isempty(nsamples)
    nsamples = floor(.05*size(D.dat,2));
  end
  
  if ~exist('snp_margin','var') || isempty(snp_margin)
    snp_margin = 500;
  end
  
  for k=1:2
    for i=1:length(regs{k})
      disp([k i]);
      h = figure();
      if show_plots
        set(h,'Visible','on');
      else
        set(h,'Visible','off');
      end
      fname = [];
      reg = regs{k}(i);
      in_chr = find(D.chrn == reg.chrn);
      chr_zero = min(in_chr);
      chr_max = max(in_chr);
      new_ads = zeros(1,length(ads{k}));
      if ~isfield(regs{k}(i),'ads')
        cur_ads = ads{k}(in_chr);
      else
        cur_ads = regs{k}(i).ads;
      end
      if k==2 && ~isempty(gene_gistic_rg)
        cur_ads = smooth_gene_scores(D,gene_gistic_rg,cur_ads,reg.chrn);
      end
      new_ads(in_chr) = cur_ads;
      range_st = max(reg.peak_wide_st-snp_margin,chr_zero);
      range_en = min(reg.peak_wide_en+snp_margin,chr_max);
      range = range_st:range_en;
      if k==1
        max_in_range = max(D.dat(range,:));
        [S si] = sort(max_in_range,'descend');
      else
        min_in_range = min(D.dat(range,:));
        [S si] = sort(min_in_range,'ascend');
      end
      
      subplot(2,4,[1:3 5:7]);
      imagesc(D.dat(range,si(1:nsamples)),[-1.5 1.5]); bluepink; 
      hold on
      wide_st = find(range == reg.peak_wide_st); wide_en = find(range == reg.peak_wide_en);
      line([1 200],[wide_st wide_st],'Color','green','LineStyle','-');
      line([1 200],[wide_en wide_en],'Color','green','LineStyle','-');
      title(['Data view for top ' num2str(nsamples) ' samples'])
      xlabel('Sample number')
      ylabel('Genomic position')
      text(-15,-35,['Chrn' num2str(reg.chrn)],'FontSize',18)
      YTickRange = 1:round(length(range)/9):length(range); 
      set(gca,'YTickLabel',{round(D.pos(range(YTickRange))/100000)*100000})
      set(gca,'YTick',YTickRange)
      
      subplot(2,4,[4 8]);
      plot(ads{k}(range),range,'Color','blue')
      hold on
      plot(new_ads(range),range,'Color','red')
      hold on
      line([0 2*max(ads{k}(range))],[reg.peak_wide_st reg.peak_wide_st],'Color','green','LineStyle','--');
      line([0 2*max(ads{k}(range))],[reg.peak_wide_en reg.peak_wide_en],'Color', ...
           'green','LineStyle','--');
      lim_s = min([ads{k}(range); new_ads(range)']);
      lim_e = max([ads{k}(range); new_ads(range)']);
      xlim([lim_s lim_e]);
      ylim([range(1) range(end)])
      set(gca,'YDir','reverse')
      title('SNP Score')
      set(gca,'YTick',[]);
            
      if isfield(reg,'genesymb') && ~isempty(reg.genesymb)
        unique_symbs = unique(reg.genesymb);
        fname = [unique_symbs{1} '.jpg'];
        for j=1:length(unique_symbs)
          gene_idx = strmatch(unique_symbs(j),{rg.symb});
          if length(gene_idx) > 1
            gene_idx = gene_idx(1);
          end
          if rg(gene_idx).chrn ~= reg.chrn
            disp(['Gene ' unique_symbs{j} ' not on same chromosome as current reg.  Skipping...']);
          else
            gene_st = rg(gene_idx).start; gene_en = rg(gene_idx).end;
            gene_snps = find_snps(D,reg.chrn,gene_st,gene_en,1);
            first_snp = min(gene_snps);
            range_idx = find(range == first_snp);
            if isempty(range_idx)
              disp(['Gene ' unique_symbs{j} ' not in range of current reg.  Skipping...']);
            else
              xl = xlim;
              text(xl(2)+.1*(xl(2)-xl(1)),first_snp,unique_symbs{j}, ...
                   'Color','red');
            end
          end
        end
      end
    
    if isempty(fname)  
      if k==1
        fname = ['Amp_reg_' num2str(i) '.jpg'];
      else
        fname = ['Del_reg_' num2str(i) '.jpg'];
      end
    end
    
    if show_plots
      keyboard
    end

    saveas(h,fname,'jpg');
    
    end
  end
