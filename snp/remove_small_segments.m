function L=remove_small_segments(L,sz,use_mb)
warning('using the smoothed data');
dat=L.smooth;
L=add_chrn(L);
for i=1:size(L.smooth,2)
  for j=1:23
    snps_in_chr=find(L.chrn==j);
    rl=runlength(L.smooth(snps_in_chr,i)');
    s212=find(rl(:,3)==1);
    if ~isempty(s212)
      for k=1:length(s212)
        st=snps_in_chr(rl(s212(k),1));
        en=snps_in_chr(rl(s212(k),2));
        if use_mb
          if (L.pos(en)-L.pos(st))<=st
            dat(st:en,i)=0;
            disp(['sample ' num2str(i) ': removing chr' num2str(j) ...
                  ': snp' num2str(st) '-snp' num2str(en)]);
          end
        else
          if (en-st+1)<=st
            dat(st:en,i)=0;
            disp(['sample ' num2str(i) ': removing chr' num2str(j) ...
                  ': snp' num2str(st) '-snp' num2str(en)]);
          end
        end
      end
    end
  end
end
