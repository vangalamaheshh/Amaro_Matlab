function X = dRanger_annotate_sites(X,P,R)
%
% given in X:
%   chr1, pos1, str1 (str=strand)
%   chr2, pos2, str2
%
% adds the following fields to X:
%   for end1 and end2 separately:
%      gene1:          name of gene (or closest gene, for IGR)
%      site1:          old-style text version "IGR: 2Mb from EGFR"
%    
%   for the rearrangement as a whole:
%      fusion:         does this rearrangement join transcribed sequenced to something else interesting?
%
% Mike Lawrence 2010-04-05

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'build_refseq',P.build);
P = impose_default_value(P,'impute_promoters',true);
P = impose_default_value(P,'imputed_promoter_size',3000);

flds = {'chr1','pos1','str1','chr2','pos2','str2'};
require_fields(X,flds);
nx = slength(X);
if ~isnumeric(X.chr1), X.chr1 = convert_chr(X.chr1);end
if ~isnumeric(X.pos1), X.pos1 = str2double(X.pos1);end
if ~isnumeric(X.str1)
  if any(strcmp(X.str1,'(+)'))
    X.str1 = listmap(X.str1,{'(+)','(-)'})-1;
  else
    X.str1 = str2double(X.str1);
  end
end
if ~isnumeric(X.chr2), X.chr2 = convert_chr(X.chr2);end
if ~isnumeric(X.pos2), X.pos2 = str2double(X.pos2);end
if ~isnumeric(X.str2)
  if any(strcmp(X.str2,'(+)'))
    X.str2 = listmap(X.str2,{'(+)','(-)'})-1;
  else
    X.str2 = str2double(X.str2);
  end
end

if ~exist('R','var'), R = load_refseq(P.build_refseq); end
if iscell(R.chr), R.chr = convert_chr(R.chr); end

R.gene_start = R.tx_start;
R.gene_end = R.tx_end;
if P.impute_promoters
  idx = find(strcmp('+',R.strand));
  R.gene_start(idx) = R.gene_start(idx) - P.imputed_promoter_size;
  idx = find(strcmp('-',R.strand));
  R.gene_end(idx) = R.gene_end(idx) + P.imputed_promoter_size;
end

% process each rearrangement

for x=1:nx
  if ~mod(x,1000), fprintf('%d/%d ',x,nx); end

  % first: identify each end

  intronnum = nan; intronframe = nan;

  for ee=1:2
    if ee==1
      chr = X.chr1(x);
      pos = X.pos1(x);
      str = X.str1(x);
    else
      chr = X.chr2(x);
      pos = X.pos2(x);
      str = X.str2(x);
    end
    thischr = find(R.chr==chr);
    overlaps = thischr(R.gene_start(thischr)<=pos & R.gene_end(thischr)>=pos);
    if isempty(overlaps)
      % it's intergenic
      dist_before = abs(R.tx_start(thischr)-pos);
      dist_after = abs(pos-R.tx_end(thischr));
      for y = 1000.^(1:0.3:3)
        idx = find(dist_before<=y);
        if ~isempty(idx)
          [tmp i] = min(dist_before(idx));
          i = idx(i);
          name = R.gene{thischr(i)};
          strand = R.strand{thischr(i)};
          zone = 0;
          a = ['IGR: ' bp2str(dist_before(i)) ' before ' name '(' strand ')'];
          break;
        end
        idx = find(dist_after<=y);
        if ~isempty(idx)
          [tmp i] = min(dist_after(idx));
          i = idx(i);
          name = R.gene{thischr(i)};
          strand = R.strand{thischr(i)};
          zone = 0;
          a = ['IGR: ' bp2str(dist_after(i)) ' after ' name '(' strand ')'];
          break;
        end
      end
    else
      % it's within a transcript: determine promoter/UTR/intron/exon
      no = length(overlaps);
      c = nan(no,1); % zone: 1=exon, 2=intron, 3=3'-UTR, 4=5'-UTR, 5=promoter
      d = nan(no,1); e = nan(no,1);  % for exons: which one, and how far
      d1 = nan(no,1); d2 = nan(no,1); e1 = nan(no,1); e2 = nan(no,1);  % for introns: between which exons and how far?
      f = nan(no,1); % for introns: how many bases in the partially completed codon?
      for j=1:no, i=overlaps(j);
        % in promoter?
        if pos<R.tx_start(i), c(j)=5; d(j) = R.tx_start(i) - pos; continue; end
        if pos>R.tx_end(i), c(j)=5; d(j) = pos - R.tx_end(i); continue; end
        % in UTR?
        if strcmp(R.strand{i},'+')
          if R.code_start(i)>pos, c(j)=4; d(j) = R.code_start(i)-pos; continue; end
          if R.code_end(i)<pos, c(j)=3; d(j) = pos-R.code_end(i); continue; end
        else % (-)
          if R.code_start(i)>pos, c(j)=3; d(j) = R.code_start(i)-pos; continue; end
          if R.code_end(i)<pos, c(j)=4; d(j) = pos-R.code_end(i); continue; end
        end
        % in exon(s)?
        in_e = find(R.exon_starts{i}<=pos & R.exon_ends{i}>=pos);
        if ~isempty(in_e), c(j)=1; e(j)=in_e; continue; end
        % otherwise: in intron
        c(j)=2;
        for k=1:R.n_exons(i)-1
          if R.exon_ends{i}(k)<pos & R.exon_starts{i}(k+1)>pos
            if strcmp(R.strand{i},'+'), f(j) = R.exon_frames{i}(k+1);
            else f(j) = R.exon_frames{i}(k);
            end
            e1(j) = k; e2(j) = k+1;
            d1(j) = pos-R.exon_ends{i}(k); d2(j) = R.exon_starts{i}(k+1)-pos;
            d(j) = min(d1(j),d2(j));           
            break;
          end         
        end
      end

      % find the transcript in the highest-priority class
      zone = -1;
      for cidx=1:5
        idx = find(c==cidx);
        if isempty(idx), continue; end
        if cidx>1
          [tmp k] = min(d(idx));
          j = idx(k);
        else
          j = idx(1);
        end
        i = overlaps(j);
        name = R.gene{i};
        strand = R.strand{i};
        if strcmp(strand,'-')
          e(j)=R.n_exons(i)-e(j)+1;
          e1(j)=R.n_exons(i)-e1(j)+1;
          e2(j)=R.n_exons(i)-e2(j)+1;
        end
        zone = cidx;
        break;
      end
      switch(zone)
       case 1  % exon
         a = sprintf('Exon %d of %s(%s)',e(j),name,strand);
       case 2  % intron
         if strcmp(strand,'+')
           if d1(j) < d2(j)
             a = sprintf('Intron of %s(%s): %s after exon %d',name,strand,bp2str(d1(j)),e1(j));
           else
             a = sprintf('Intron of %s(%s): %s before exon %d',name,strand,bp2str(d2(j)),e2(j));
           end
         else % (-)
           if d1(j) < d2(j)
             a = sprintf('Intron of %s(%s): %s before exon %d',name,strand,bp2str(d1(j)),e1(j));
           else
             a = sprintf('Intron of %s(%s): %s after exon %d',name,strand,bp2str(d2(j)),e2(j));
           end
         end
         intronnum = e1(j);
         intronframe = f(j);
       case 3  % 3'-UTR
         a = sprintf('3''-UTR of %s(%s): %s after coding stop',name,strand,bp2str(d(j)));
       case 4  % 5'-UTR
         a = sprintf('5''-UTR of %s(%s): %s before coding start',name,strand,bp2str(d(j)));
       case 5  % promoter
         a = sprintf('Promoter of %s(%s): %s from tx start',name,strand,bp2str(d(j)));
       case -1
        a = sfprintf('unexpected error: zone = -1');
      end
    end

    if ee==1
      gene1 = name;
      X.gene1{x,1} = gene1;
      X.site1{x,1} = a;      
      zone1 = zone;
      strandgene1 = strand;
      str1 = str;
      intronnum1 = intronnum;
      intronframe1 = intronframe;
    else
      gene2 = name;
      X.gene2{x,1} = gene2;
      X.site2{x,1} = a;
      zone2 = zone;
      strandgene2 = strand;
      str2 = str;
      intronnum2 = intronnum;
      intronframe2 = intronframe;
    end

  end   % next end of this rearrangement

  % does this rearrangement generate an interesting fusion product?
  %    zone: 0=IGR 1=exon, 2=intron, 3=3'-UTR, 4=5'-UTR, 5=promoter
  
  strandmatch1 = ((strcmp(strandgene1,'+') & str1==0) | (strcmp(strandgene1,'-') & str1==1));
  strandmatch2 = ((strcmp(strandgene2,'+') & str2==0) | (strcmp(strandgene2,'-') & str2==1));
  txactive1 = (zone1>0 && zone1~=3);
  txactive2 = (zone2>0 && zone2~=3);
  
  txt = '-';
  if (txactive1 && strandmatch1) && (txactive2 && strandmatch2)
    txt = 'Antisense fusion';
  elseif (txactive1 && strandmatch1) || (txactive2 && strandmatch2)
    if txactive1 && txactive2
      if strcmp(gene1,gene2) % within-gene event      
        if str1==0 && str2==1
          class = 'Deletion';
        elseif str1==1 && str2==0
          class = 'Duplication';
        else
          class = 'Inversion';
        end
        if zone1==2 && zone2==2         
          if intronnum1==intronnum2
            txt = [class ' within intron'];
          else
            numexons = abs(intronnum2-intronnum1);
            txt = [class ' of ' num2str(numexons) ' exon'];
            if numexons>1, txt = [txt 's']; end
            if ~strcmpi(class,'Inversion')
              if intronframe1==intronframe2
                txt = [txt ': in frame'];
              else
                txt = [txt ': out of frame'];
              end
            end
          end
        else
          txt = [class ' within transcript'];
          if zone1==1 || zone2==1
            txt = [txt ': mid-exon'];
          end
        end        
      else % inter-gene event
        if strandmatch1, fusname = [gene1 '-' gene2]; else fusname = [gene2 '-' gene1]; end
        if zone1==2 && zone2==2
          if intronframe1==intronframe2
            txt = 'Protein fusion: in frame';
          else
            txt = 'Protein fusion: out of frame';
          end
        elseif zone1==1 || zone2==1
          txt = 'Protein fusion: mid-exon';
        else
          txt = 'Transcript fusion';
        end
        txt = [txt ' (' fusname ')'];
      end  
    end
  end
  X.fusion{x,1} = txt;
    
end   % next rearrangement

fprintf('\nDone\n');

