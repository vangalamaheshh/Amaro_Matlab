function regs = add_genes(D,regs,cyto,rg,partial_hits)
%ADD_GENES Add gene lists to a lists of amplification and deletion peaks
%
%   REGS = add_genes(D,REGS,CYTO,RG,PARTIAL_HITS)
%
%   D is a copy number data structure
%   REGS is a two-element amp/del cell array containing arrays of structs
% for amplification and deletion peaks. ADD_GENES adds a 'genes' field to
% each struct array containing the genes associated with the peak.
%   CYTO is the struct array containing cytoband information
%   RG is the reference genome struct array
%   PARTIAL_HITS is an optional logical parameter that controls whether
% genes partially in the peak are included (1) or not (0). The default is
% 1. If PARTIAL_HITS is a vector, then it is treated as an amp/del pair of
% values for treating amplifications deletions with different stringencies.
%

  % process optional, polymorphic partial_hits parameter
  if ~exist('partial_hits','var') || isempty(partial_hits)
    partial_hits = 1;
  end
  if isscalar(partial_hits)
    partial_hits = [partial_hits partial_hits];
  end
  
  % amps, then dels
  for k=1:2
    % loop over peaks
    for i=1:length(regs{k})
      verbose([num2str(k) ' ' num2str(i)],30);
      cur_reg = regs{k}(i);
      [~,~,enl_region_to_next_snp] = genomic_location(D, ...
                    {[cur_reg.peak_wide_st,cur_reg.peak_wide_en]},cyto,1);
      genes_in_region = genes_at(rg,cur_reg.chrn, ...
                                 enl_region_to_next_snp{1}(1),...
                                 enl_region_to_next_snp{1}(2),...
                                 1,partial_hits(k));
      if ~isempty(genes_in_region)
        lid_in_region = cat(1,rg(genes_in_region).locus_id);
        [~,un_lid_reg_i1,~] = unique(lid_in_region);
        genes_in_region = genes_in_region(un_lid_reg_i1);
      end
      
      regs{k}(i).genes = {rg(genes_in_region).symb};
    end
  end
