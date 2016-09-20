function [as,names] = annotate_site(D,chrs,sts,ens)
%
% OBSOLETE:  replaced by dRanger_annotate_sites
%
% for a given chr:st-en region, looks up region in passed database D
% (as loaded by load_genedb)
% and returns a 'genomic_position' annotation for end1 and for end2.
%
% example annotations (listed in order of priority):
%
% 1.  "Exon 3 of EGFR"
% 2.  "3'-UTR of EGFR: 10bp from last exon"
% 3.  "5'-UTR of EGFR: 20bp from first exon"
% 4.  "Intron of EGFR: 6bp from exon 1" (gives distance to nearest exon)
% 5.  "IGR: 2Mb from EGFR" (gives distance to closest gene)
%
% Also returns name of closest gene as "name"
%
% Mike Lawrence 2009-02-10

% STEPS:
%   (1) compute distance from each transcript in database (negative numbers indicate overlap)
%   (2) if none overlap (i.e. class 5, IGR)
%       --> pick the closest one and use that
%   (3) otherwise:
%       --> find out what class each overlap is in (1/2/3/4),
%           and use the one in the highest-priority class.

if iscell(D.chrom)
  fprintf('Converting D.chrom to numeric...\n');
  D.chrom = convert_chr(D.chrom);
end
Dfull = D;

if ~exist('ens','var'), ens=sts; end

as = cell(length(chrs),1);
names = cell(length(chrs),1);
for iii=1:length(chrs)
chr=chrs(iii); st=sts(iii); en=ens(iii);

D = reorder_struct(Dfull,Dfull.chrom==chr);

overlaps = find(D.chrom==chr & D.txStart<en & D.txEnd>st);

if isempty(overlaps)
% it's intergenic
  dist = min(abs(D.txStart-en),abs(st-D.txEnd));
  for x = 1000.^(1:0.3:3)
    idx = find(dist<=x);
    if ~isempty(idx)
      idx2 = idx(find(strcmp(D.db(idx),'RefSeq')));    % prefer RefSeq because they have real gene names
      if isempty(idx2), idx2 = idx; end
      [tmp i] = min(dist(idx2));
      i = idx2(i);
      if strcmp(D.db{i},'RefSeq'), name = D.name2{i};
      else name = D.name{i}; end;
      strand = D.strand{i};
      a = ['IGR: ' bp2str(dist(i)) ' from ' name '(' strand ')'];
      break;
     end
  end
else

% it's within a transcript

no = length(overlaps);
c = nan(no,1); d = nan(no,1); e1 = nan(no,1); e2 = nan(no,1);
for j=1:no,i=overlaps(j);
  % in exon(s)?
  in_e = find(D.exonStarts{i}<=en & D.exonEnds{i}>=st);
  if ~isempty(in_e), c(j)=1; e1(j)=min(in_e); e2(j)=max(in_e); continue; end
  % in UTR?
  if strcmp(D.strand{i},'+')
    if D.cdsStart(i)>en, c(j)=3; d(j) = D.cdsStart(i)-en; continue; end
    if D.cdsEnd(i)<st, c(j)=2; d(j) = st-D.cdsEnd(i); continue; end
  else % (-)
    if D.cdsStart(i)>en, c(j)=2; d(j) = D.cdsStart(i)-en; continue; end
    if D.cdsEnd(i)<st, c(j)=3; d(j) = st-D.cdsEnd(i); continue; end
  end
  % otherwise: in intron
  c(j)=4;
  for e=1:D.exonCount(i)-1
    if D.exonEnds{i}(e)<st & D.exonStarts{i}(e+1)>en
      if st-D.exonEnds{i}(e) < D.exonStarts{i}(e+1)-en % closer to left exon
        e1(j)=e; d(j)=st-D.exonEnds{i}(e);
      else e1(j)=e+1; d(j)=D.exonStarts{i}(e+1)-en; end
    end
  end
  if isnan(d(j))
    fprintf('Inconsistent behavior #1 in annotate_site\n');
    disp(c(j));
  end
end

% find the transcript in the highest-priority class

for cidx=1:4
  idx = find(c==cidx);
  if isempty(idx), continue; end
  idx2 = idx(find(strcmp(D.db(overlaps(idx)),'RefSeq')));
  if isempty(idx2), idx2 = idx; end
  if cidx>1
    [tmp j] = min(d(idx2));
  else
    j=1;
  end
  i = overlaps(idx2(j));
  if strcmp(D.db{i},'RefSeq'), name = D.name2{i};
  else name = D.name{i}; end;
  strand = D.strand{i};
  if strcmp(strand,'-')
    e1(j)=D.exonCount(i)-e1(j)+1;
    e2(j)=D.exonCount(i)-e2(j)+1;
  end
  switch(cidx)
    case 1  % exon
      if e2(j)>e1(j), a = sprintf('Exons %d-%d of %s(%s)',e1(j),e2(j),name,strand);
      else a = sprintf('Exon %d of %s(%s)',e1(j),name,strand); end
    case 2  % 3'-UTR
      a = sprintf('3''-UTR of %s(%s): %s from last exon',name,strand,bp2str(d(j)));
    case 3  % 5'-UTR
      a = sprintf('5''-UTR of %s(%s): %s from first exon',name,strand,bp2str(d(j)));
    case 4  % intron
      a = sprintf('Intron of %s(%s): %s from exon %d',name,strand,bp2str(d(j)),e1(j));
  end
  break;
end

end

as{iii} = a;
names{iii} = name;
end

if length(chrs)==1
  as = a;
  names = names;
end
