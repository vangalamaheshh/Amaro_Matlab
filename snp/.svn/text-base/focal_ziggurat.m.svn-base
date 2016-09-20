function [QA QD QAOD QDOA] = focal_ziggurat(D,cyto,broad_levels,centromere_snps)

%% Focal_ziggurat performs ziggurat peel-off.  It has been modified
%from the original ziggurat.m function to  allow ziggurat to be performed with respect
%to an underlying broad_level(s) on each chromosome.  These
%broad_levels are represented in the broad_levels matrix, with
%corresponding 'breakpoints' in the centromere_snps matrix.
  
%broad_levels = 2n x j matrix where n=# chromosomes, j = # samples
%broad_levels(2i-1,j) = broad_level on 'p-arm' of i-th chr in j-th sample
%broad_levels(2i,j) = broad_level on 'q-arm' of i-th chr in j-th sample  
  
%centromere_snps = n x j matrix where n=# chromosomes, j= # samples
%centromere_snps(i,j) = location of broad breakpoint on i-th chr in
%j-th sample [Note that this may, but is not generally, the actual
%centromere on the i-th chromosome]
  
%For chromosomes with no break, centromere_snps(i,j) = NaN
%In these cases, the broad_level for that chromosome is located in broad_levels(2i,j)
    
 
  changeChrn = diff(D.chrn);
  chrnEnd = find(changeChrn == 1);
  chrnEnd(end+1) = size(D.chrn,1);
  chrnStart(1)=1;
  chrnStart(2:length(chrnEnd)) = chrnEnd(1:end-1)+1;
  
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
  
  armlengths = [chrarms{1}.length];
  armstart = [chrarms{1}.start];
  
  %% find segments
  diffDat = diff(D.dat);
  [breakpt_ids sample_ids] = find(diffDat ~= 0);
    
  for j=1:length(unique(sample_ids))
    disp(j)
    with_cent = find(isnan(centromere_snps(:,j)) == 0);
    no_cent = setdiff(1:max(D.chrn),with_cent)';
    represented_arms = sort([2*with_cent-1; 2*with_cent; 2*no_cent]);
    bpt = breakpt_ids(find(sample_ids == j));
    bpt = union(bpt,chrnEnd);
    if ~isempty(with_cent)
      bpt = union(bpt,centromere_snps(with_cent,j));
    end
    B = zeros(length(bpt),5);
    B(:,5) = repmat(j,length(bpt),1);
    B(1,2) = 1;
    B(2:end,2) = bpt(1:end-1)+1;
    B(:,3)=bpt;
    B(:,1)=D.chrn(B(:,2));
    B(:,4)=D.dat(B(:,2),j);
    B(:,4)=(2.^(B(:,4)+1))-2;
    sample_broad_levels = (2.^(broad_levels(represented_arms,j)+1))-2;
    [Qamp{j} Qdel{j} Qampondel{j} Qdelonamp{j}] = focal_ziggurat_deconstruction(B,sample_broad_levels,with_cent,centromere_snps(with_cent,j));
    
  end

  QA = cat(1,Qamp{:}); 
  QD = cat(1,Qdel{:});
  QAOD = cat(1,Qampondel{:});
  QDOA = cat(1,Qdelonamp{:});
