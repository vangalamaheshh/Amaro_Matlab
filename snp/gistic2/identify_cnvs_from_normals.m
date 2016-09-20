function identify_cnvs_from_normals(outfile,D,pop_cnv_file,normal_al,amplitude_p,min_fract,remove_X)

          
     %% Get normal arrays
     verbose('Subselecting normals...',30);
     [M m1 m2] = match_string_sets_hash(D.sdesc,normal_al); 
     verbose(['Found ' num2str(length(m2)) ' of ' num2str(length(normal_al)) ' samples'],30);
     Dn = reorder_D_cols(D,m1);
     
     % Remove X chromosome
     if remove_X
       verbose('Removing X chromosome...',30);
       Dnx = reorder_D_rows(Dn,find(Dn.chrn <= 22));
     else
       Dnx = Dn;
     end
     
     %% Find cut-offs for normals
     verbose('Finding cutoffs...',30);
     hx = -2:0.01:2;
     hc = sum(histc(Dnx.dat,hx),2);
     
     ps = hc/sum(hc);
     [m mid] = max(ps);
     win_size = 0;
     comb_p = sum(ps(mid-win_size:mid+win_size));
     while comb_p < (1-amplitude_p)
       win_size = win_size+1;
       comb_p = sum(ps(mid-win_size:mid+win_size));
     end
                 
     lcutoff = hx(mid-win_size+1)
     ucutoff = hx(mid+win_size+1)

     %% Find snps altered in greater than min_fract% of normals
     verbose('Finding cnvs in normals...',30);
     qq = (Dnx.dat >= ucutoff) + (Dnx.dat <= lcutoff);
     nl_cnvs = find(sum(qq,2) >= (min_fract*size(Dnx.dat,2)));
     
     %% Find population cnvs
     verbose('Finding population cnvs...',30);
     [Drm,nsnps,rm_snps] = remove_cnv(Dnx,pop_cnv_file);
     pop_cnvs = find(rm_snps);
     
     % Merge into one big list of cnvs
     verbose('Merging lists...',30);
     cnvs = union(nl_cnvs,pop_cnvs);
     
     vv = zeros(1,size(qq,1));
     vv(cnvs) = 1;
     rl = runlength(vv);
     
     cnv_idx = find(rl(:,3) == 1);
     cnv_snps = rl(cnv_idx,1:2);
     
     %% Fix those that cross chromosome boundary
     cross_chr = find(Dnx.chrn(cnv_snps(:,2))-Dnx.chrn(cnv_snps(:, ...
                                                       1)));
     
     for j=1:length(cross_chr)
       first_chr = Dnx.chrn(cnv_snps(cross_chr(j),1));
       second_chr = Dnx.chrn(cnv_snps(cross_chr(j),2));
       first_chr_end = max(find(Dnx.chrn == first_chr));
       second_chr_st = min(find(Dnx.chrn == second_chr));
       tmp = cnv_snps(cross_chr(j),2);
       cnv_snps(cross_chr(j),2) = first_chr_end;
       cnv_snps(size(cnv_snps,1)+1,:) = [second_chr_st tmp];
     end
     
     %% sort
     [S si] = sort(cnv_snps(:,1));
     cnv_snps = cnv_snps(si,:);
          
     %% Now write out cnv file
     verbose('Writing cnv file...',30);
     f = fopen(outfile,'w')

     fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\n','ID','Chromosome','Start', ...
        'End','Flanking start','Flanking end');

     for j=1:size(cnv_snps,1)
       modi(j,100)
       ids = ['CNP-' num2str(j)];
       chrn = Dnx.chrn(cnv_snps(j,1));
       st = Dnx.pos(cnv_snps(j,1));
       en = Dnx.pos(cnv_snps(j,2));
       %% Check if at start of chromosome (if so, flanking pos = 1)
       if min(find(Dnx.chrn == chrn)) == cnv_snps(j,1)
         fl_st = 1;
       else
         fl_st = Dnx.pos(cnv_snps(j,1)-1);
       end
       %% Check if at end of chromosome (if so, flanking pos = end of chrn)
       if max(find(Dnx.chrn == chrn)) == cnv_snps(j,2)
         fl_en = en+1;
       else
         fl_en = Dnx.pos(cnv_snps(j,2)+1);
       end
       fprintf(f,'%s\t%d\t%d\t%d\t%d\t%d\n',ids,chrn,st,en,fl_st,fl_en);
     end
     
     fclose(f)
     




     
     
       
       
    