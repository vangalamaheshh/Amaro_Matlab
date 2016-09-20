function new_dat = subtract_broad_level_from_chrarm(CL,cyto,broad_levels)

  %% Broad levels is a 44 x n matrix, where n is the number of samples
  %% Each element represents the broad level for that chromosome arm in
  %% that sample  
  
  new_dat = CL.dat;
  
  armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9]+[p-q]+','match'),'UniformOutput',false));
  armnames = armnames(2:end);

  chrarms = {};
  for i=1:length(armnames)
    idx = strmatch(armnames(i),{cyto.name});
    band.name = char(armnames(i));
    band.start = cyto(idx(1)).start+1;
    band.end = cyto(idx(end)).end;
    band.chrn=cyto(idx(end)).chrn;
    band.length = band.end-band.start+1;
    chrarms{1}(i) = band;
  end

  chr=cat(1,chrarms{1}.chrn);
  [spos,sposi]=sort(chr);
  chrarms{1} = chrarms{1}(sposi);
  
  for i=1:22
    centromeres(i) = chrarms{1}(2*i-1).end;
  end

  disp('Subtracting broad level from each arm on chromosome:')
  for ch=1:22
    disp(ch)
    ch_snps = find(CL.chrn == ch);
    chr_snp_start = min(ch_snps);
    chr_snp_cent = ch_snps(max(find(CL.pos(ch_snps)< centromeres(ch))));
    chr_snp_end = max(ch_snps);
    if ~isempty(chr_snp_cent)
      new_dat(chr_snp_start:chr_snp_cent,:) = new_dat(chr_snp_start:chr_snp_cent,:)-repmat(broad_levels(2*ch-1,:),length(chr_snp_start:chr_snp_cent),1);
      new_dat(chr_snp_cent+1:chr_snp_end,:) = new_dat(chr_snp_cent+1:chr_snp_end,:)-repmat(broad_levels(2*ch,:),length(chr_snp_cent+1:chr_snp_end),1);
    else
      new_dat(chr_snp_start:chr_snp_end,:) = new_dat(chr_snp_start:chr_snp_end,:)-repmat(broad_levels(2*ch,:),length(chr_snp_start:chr_snp_end),1);
    end
  end
  
      
