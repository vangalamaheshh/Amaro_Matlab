function M = get_allele_fracs(M,P)
% get_allele_fracs(M,P)
%   demand_fields(M,{'tbam','chr','start','end','ref_allele','newbase'});

if ~exist('P','var'), P=[]; end

if ~isfield(P,'build') && ~isfield(M,'build')
  fprintf('Assuming build hg18\n');
  P.build = 'hg18';
end

P = impose_default_value(P,'minqual',20);
P.include_unmapped_reads = false;

demand_fields(M,{'tbam','chr','start','end','ref_allele','newbase'});
% tbam = location of tumor bam
M = make_numeric(M,{'start','end'});
M.chr = convert_chr(M.chr);

demand_files(M.tbam);
% make sure BAM indeces are also available
M.tbai = regexprep(M.tbam,'\.bam$','.bai');
ok = demand_file(M.tbai);
if any(~ok)
  tmp = M.tbai(~ok);
  tmp = regexprep(tmp,'\.bai$','.bam.bai');
  demand_files(tmp);
  fprintf('But *.bam.bai files ARE available in all of the above cases.\n');
end

refbase('ACGT')=1:4;
altbase('ACGT-')=[65:68 -100];

M.tot_count = nan(slength(M),1);
M.alt_count = nan(slength(M),1);
M.ref_count = nan(slength(M),1);
M.allele_frac = nan(slength(M),1);

for i=1:slength(M)
  fprintf('%d/%d   ',i,slength(M));
  if isfield(M,'patient_name'), fprintf('%s   ',M.patient_name{i});
  elseif isfield(M,'Tumor_Sample_Barcode'), fprintf('%s   ', M.Tumor_Sample_Barcode{i});
  elseif isfield(M,'patient') && ischar(M.patient{i}), fprintf('%s   ', M.patient{i});
  end
  if isfield(M,'gene_name'), fprintf('%s   ',M.gene_name{i});
  elseif isfield(M,'Hugo_Symbol'), fprintf('%s   ', M.Hugo_Symbol{i});
  elseif isfield(M,'gene') && ischar(M.gene{i}), fprintf('%s   ', M.gene{i});
  end
  fprintf('chr%d:%d-%d   %s -> %s\n', M.chr(i),M.start(i),M.end(i),M.ref_allele{i},M.newbase{i});
  try
    if isfield(M,'build')
      if isnumeric(M.build), P.build = M.build(i);
      else P.build = M.build{i};
      end
    end
    [R B S] = pull_from_bam(M.tbam{i},M.chr(i),M.start(i),M.end(i),P);
    if strcmp(M.ref_allele{i},'-')
      fprintf('    Insertion: can''t tell allele frac from basecall data\n');
      readlen = R(:,5)-R(:,4)+1;
      lenmode = getBAMFileReadLength(M.tbam{i});
      inslen = length(M.newbase{i});
      fprintf('    Approximating by counting number of reads %d bp shorter than readlength mode %d\n',...
            inslen,lenmode);
      refcount = sum(readlen==lenmode);
      altcount = sum(readlen==lenmode-inslen);
      totcount = length(readlen);
    else
      B = B(B(:,4)==M.start(i),:);
      B = B(B(:,2)>=P.minqual | B(:,2)==-100,:);
      refi = refbase(M.ref_allele{i}(1));
      alti = altbase(M.newbase{i}(1));
      refcount = sum(B(:,1)==refi);
      altcount = sum(B(:,1)==alti);
      totcount = size(B,1);
    end
    M.ref_count(i) = refcount;
    M.alt_count(i) = altcount;
    M.tot_count(i) = totcount;
    M.allele_frac(i) = (altcount / totcount);
    fprintf('    %d/%d (%0.2f) reads have the alternate allele\n', ...
              altcount, totcount, altcount/totcount);
  catch me
    fprintf('    Exception!\n');
    disp(me);
    disp(me.message);
  end
end


