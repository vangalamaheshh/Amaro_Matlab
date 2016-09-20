function M = survey_refseq_transcripts(build)
% survey_refseq_transcripts(build)
%
% go through transcripts and find any with problems (specifically, start/end don't match frame)
%
% Mike Lawrence 2010-10-01

if nargout==1
  genmuts = true;
%  genmuts_in_which_exon = +2;
%  genmuts_in_which_exon = +1;
  genmuts_in_which_exon = 0;
  rand('twister',456);
  M = []; M.build = []; M.chr = []; M.start = [];
elseif nargout==0
  genmuts = false;
else
  error('too many outputs expected');
end

R = load_refseq(build);

[g gi gj] = unique(R.gene);
for gno=1:length(g)
  idx = find(gj==gno);
  for j=1:length(idx), i=idx(j);
    plusstrand = strcmp(R.strand{i},'+');
    if plusstrand, forfrom=1; forstep=+1; forto=R.n_exons(i);
    else forfrom=R.n_exons(i); forstep=-1; forto=1;
    end   
    orflen = 0;
    frameshifts = 0;
    frameshift_exon = nan;
    last_e = nan;
    for e=forfrom:forstep:forto
      st = R.exon_starts{i}(e);
      en = R.exon_ends{i}(e);
      fr = R.exon_frames{i}(e);
      if orflen>0 && fr~=-1 && mod(orflen,3)~=fr
        % recover from known programmed frameshift somewhere in the preceding exon
        frameshifts=frameshifts+1;
        frameshift_exon = last_e;
        fprintf('%s\t%s\t%d\t%s\t%d\t%d\n',R.gene{i},R.transcript{i},length(idx),R.chr{i},st,en);
        pfs = mod(fr-orflen,3);
        if pfs==1, orflen=orflen+1;
        else orflen=orflen+2;
        end
      end
      if st<=R.code_end(i) && en>=R.code_start(i)
        if R.code_start(i)>st, st = R.code_start(i); end
        if R.code_end(i)<en, en = R.code_end(i); end
        orflen = orflen + (en-st+1);
      end
      last_e = e;
    end  % next exon
    if genmuts==true && frameshifts>0
      % generate a mutation in the middle of an exon near the frameshift
      if plusstrand, e = frameshift_exon + genmuts_in_which_exon;
      else e = frameshift_exon - genmuts_in_which_exon;
      end
      if e>=1 && e<=R.n_exons(i)
        M.chr{end+1,1} = R.chr{i};
        M.start(end+1,1) = round((R.exon_starts{i}(e)+R.exon_ends{i}(e))/2);
      end
    end

  end % next transcript
end % next gene

if genmuts
  M.build = repmat({build},length(M.chr),1);
  M.end = M.start;
  chrno = convert_chr(M.chr);
  [u ui] = unique([chrno M.start],'rows');
  M = reorder_struct(M,ui);
  M = reorder_struct(M,~isnan(chrno(ui)));
  nm = slength(M);
  bases = 'ACGT';
  for i=1:nm
    M.ref_allele{i,1} = upper(genome_region(M.chr(i),M.start(i),M.end(i),build));
    newbases = setdiff(bases,M.ref_allele{i,1});
    M.tum_allele1{i,1} = newbases(ceil(3*rand));
  end
  M.tum_allele2 = M.tum_allele1;
  M.tumor_barcode = repmat({'test_tumor'},nm,1);
  M.normal_barcode = repmat({'test_normal'},nm,1);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M18=survey_refseq_transcripts('hg18');
save_struct_noheader(M18,'/xchip/cga1/lawrence/annot/test_muts_hg18.maflite');
M19=survey_refseq_transcripts('hg19');
save_struct_noheader(M19,'/xchip/cga1/lawrence/annot/test_muts_hg19.maflite');
annotate_maflite('/xchip/cga1/lawrence/annot/test_muts_hg18.maflite',['/xchip/cga1/lawrence/annot/' ...
                    'test_muts_hg18.maf.annotated'],'hg18');
annotate_maflite('/xchip/cga1/lawrence/annot/test_muts_hg19.maflite',['/xchip/cga1/lawrence/annot/' ...
                    'test_muts_hg19.maf.annotated.no_frames'],'hg19');

% which_exon = 0
M19=survey_refseq_transcripts('hg19');
save_struct_noheader(M19,'/xchip/cga1/lawrence/annot/test_muts_same_exon_hg19.maflite');
annotate_maflite('/xchip/cga1/lawrence/annot/test_muts_same_exon_hg19.maflite',['/xchip/cga1/lawrence/annot/' ...
                    'test_muts_same_exon_hg19.maf.annotated'],'hg19');
