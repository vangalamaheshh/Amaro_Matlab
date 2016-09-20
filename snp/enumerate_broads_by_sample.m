function [broad_levels breakpoints] = enumerate_broads_by_sample(D,cyto,centromere_snps)
  
  for ch = 1:max(D.chrn) %% ch represents each chromosome number
    disp(['Finding broads on chromosome ' num2str(ch)])
    chr_snps = find(D.chrn == ch);
    cent_snp = centromere_snps(ch);
    for j=1:size(D.dat,2) %% j represents each sample
      disp(j)
      cur_chr = D.dat(chr_snps,j);
      bpts = find(diff(cur_chr)~=0);
      bpts = union(bpts,max(chr_snps));
      B=zeros(length(bpts),5);
      B(:,1) = repmat(ch,length(bpts),1);
      B(:,5) = repmat(j,length(bpts),1);
      B(1,2) = min(chr_snps);
      B(2:end,2) = bpts(1:end-1)+1;
      B(:,3)=bpts;
      B(:,4)=D.dat(B(:,2),j);
      B(:,4)=(2.^(B(:,4)+1))-2;
      
      
      %% case 1: no broad level
      unique_levels = union(unique(B(:,4)),0);
      with_cent = 0;
      for i=1:length(unique_levels)
        sample_broad_levels = unique_levels(i);
        [Qamp1{i} Qdel1{i} Qampondel1{i} Qdelonamp1{i}] = ...
            focal_ziggurat_deconstruction(B,sample_broad_levels, ...
                                          with_cent,[]);
        %% How to score each one?
      end
      
      %% case 2: one broad level
      for k=1:length(bpts)-1
        with_cent = ch;
        parm = [B(1,1) min(chr_snps) bpts(k) B(1,4:5)];
        fract = normalize_by_arm_length(D,parm,cyto,1,2);
        if fract < 1 %% broad is on q-arm
          idx = min(find(B(:,2)>=bpts(k)));
          unique_levels = setdiff(unique(B(idx:end,4)),0);
          si = [0 1];
        else %% Broad is on p-arm
          idx = find(B(:,3)==bpts(k));
          unique_levels = setdiff(unique(B(1:idx,4)),0);
          si = [1 0];
        end
        
        for l=1:length(unique_levels)
          sample_broad_levels = [unique_levels(l)*si(1) unique_levels(l)* ...
                              si(2)];
          [Qamp2{k,l} Qdel2{k,l} Qampondel2{k,l} Qdelonamp2{k,l}] = ...
              focal_ziggurat_deconstruction(B,sample_broad_levels,with_cent,bpts(k));
          
          %% how to score each one?
          
        end
        
      end
      
      %% case 3: two broad levels centered at centromere
      with_cent = ch;
      cent_idx = find(B(:,2)<=cent_snp & B(:,3)>=cent_snp);
      if B(cent_idx,2) == cent_snp %% seg at cent_idx starts at centromere
        plevels = setdiff(unique(B(1:cent_idx-1,4)),0);
        qlevels = setdiff(unique(B(cent_idx:end,4)),0);
        cent_idx = cent_idx-1;
      elseif B(cent_idx,3) == cent_snp %% seg at cent_idx ends at
                                       %centromere
        plevels = setdiff(unique(B(1:cent_idx,4)),0);
        qlevels = setdiff(unique(B(cent_idx+1:end,4)),0);
      else %% snp at cent_idx spans centromere
        tmp = B(cent_idx,:);
        B(cent_idx,3) = cent_snp;
        tmp(2) = cent_snp+1;
        B(cent_idx+2:size(B,1)+1,:) = B(cent_idx+1:end,:);
        B(cent_idx+1,:) = tmp;
        plevels = setdiff(unique(B(1:cent_idx,4)),0);
        qlevels = setdiff(unique(B(cent_idx+1:end,4)),0);
      end
      
      for l = 1:length(plevels)
        sample_broad_levels = [plevels(l) 0];
        [Qamp3p{l} Qdel3p{l} Qampondel3p{l} Qdelonamp3p{l}] = ...
              focal_ziggurat_deconstruction(B,sample_broad_levels, ...
                                            with_cent,cent_snp);
      end
      
      for l=1:length(qlevels)
         sample_broad_levels = [0 qlevels(l)];
         [Qamp3q{l} Qdel3q{l} Qampondel3q{l} Qdelonamp3q{l}] = ...
              focal_ziggurat_deconstruction(B,sample_broad_levels, ...
                                            with_cent,cent_snp);
      end
    
      %% How to score
    end
  end
      
