function summary = write_gene_summary(cur_rg,num_tumors,types,data_mat,...
        sig_types_per_gene,num_sig_peaks_per_gene,isamp,isubtypes,all,qv_thresh)
  
  if ~exist('qv_thresh','var') || isempty(qv_thresh)
      qv_thresh = 0.25;
  end
  in_peak = logical(data_mat(:,1))';
  num_genes = data_mat(:,5);
  qvs = data_mat(:,6);
      
  % Find all_tumors
  cur_gene = cur_rg.symb;
  
  if isamp
    desc_str = 'amplified';
    desc_str1 = 'amplification';
  else
    desc_str = 'deleted';
    desc_str1 = 'deletion';
  end
  
  if qvs(all) <= qv_thresh
    summary = [cur_gene ' is significantly focally ' desc_str ' across the entire dataset of ' num2str(num_tumors) ' tumors ']; 
    if in_peak(all)
      summary = [summary 'and is located within a focal peak region of ' desc_str1 ...
                 ' containing ' num2str(num_genes(all)-1) ' additional genes.'];
    else % There must be a significant peak on chromosome b/c qv < .25
      summary = [summary 'but is not located within a focal peak region of ' ...
                 desc_str1 '.'];
      peak_start = data_mat(all,3);
      peak_end = data_mat(all,4);
      if cur_rg.start > peak_end
        distance = double(cur_rg.start-peak_end);
      else
        distance = double(peak_start-cur_rg.end);
      end
      summary = [summary ' ' cur_gene ' is located ' ...
                 num2str(floor(distance/1e4)/100) ' Mb away from the nearest peak region of ' desc_str1 '.'];
    end
  else % qv across dataset is greater than .25
    summary = [cur_gene ' is not significantly focally ' desc_str ...
               ' across the entire dataset of ' num2str(num_tumors) ' tumors ']; 
    if in_peak(all)
      summary = [summary '.  However, it is located within a focal peak region of ' desc_str1 ...
                 ' containing ' num2str(num_genes(all)-1) ' additional genes ( '];
      summary = [summary 'this is most likely due to poor peak definition)'];
    else
      summary = [summary 'and is not located within a focal peak region of ' ...
                desc_str1 '.'];
      if ~isnan(data_mat(all,3))
        peak_start = data_mat(all,3);
        peak_end = data_mat(all,4);
        if cur_rg.start > peak_end
          distance = double(cur_rg.start-peak_end);
        else
          distance = double(peak_start-cur_rg.end);
        end
        summary = [summary ' ' cur_gene ' is located ' num2str(floor(distance/1e4)/100) ...
                   ' Mb away from the nearest peak region of ' desc_str1 '.'];
      end
    end     
  end
    
  if ~exist('isubtypes','var') || isempty(isubtypes)
    isubtypes = true(size(types));
    isubtypes(all) = false;
  end
      
  sig_types = (qvs <= qv_thresh)' & isubtypes;
%!  sig_types = find(qvs(idx) <= .25);

  
  if nnz(sig_types) >= 1
    summary = [summary '  ' cur_gene ' is significantly focally ' desc_str ...
               ' in ' num2str(nnz(sig_types)) ' of ' num2str(nnz(isubtypes)) ...
               ' independent subtypes analyzed in our dataset.'];
    also_in_peak = sig_types & in_peak;
%    also_in_peak = sig_types(find(in_peak(sig_types)));
    if nnz(also_in_peak)
      summary = [summary ' Among these, it is located within a focal peak region of ' ...
                 desc_str1 ' in ' num2str(nnz(also_in_peak)) ' subtypes.']; 
    else
      summary = [summary ' It is not within the focal peak region of ' ...
                 desc_str1 ' in any of the individual tumor types.'];
    end
    fract_sig = sum(sig_types_per_gene >= nnz(sig_types)) / ...
                                            length(sig_types_per_gene);
    fract_peaks = sum(num_sig_peaks_per_gene >= nnz(also_in_peak)) / ...
                                            length(num_sig_peaks_per_gene);
    summary = [summary ' For reference, ' sprintf('%0.3f%',100*fract_sig) ...
               '% of all genes are significantly focally ' desc_str ...
               ' in at least '  num2str(nnz(sig_types)) ' subtypes and ' ...
               sprintf('%0.3f%',100*fract_peaks) '% of all genes are present in focal ' ...
               desc_str1 ' peaks in at least ' num2str(nnz(also_in_peak)) ' subtypes.'];
  else
    summary = [summary ' ' cur_gene ' is not significantly focally ' ...
               desc_str ' in any of the ' num2str(nnz(isubtypes)) ...
               ' individual subtypes analyzed in our dataset.'];
  end
  
  %summary = [summary ' When considering the significance of any individual gene, additional factors to take into account include: the number of genes in each peak, the proximity to known cancer genes, and the significance (q-value) of the peak.'];
  
  
