function X = dRanger_check_primers(X,P)
% Mike Lawrence 2010-04-12

if ~exist('P','var'), P = []; end

P = impose_default_value(P,'SNP_mask_BAMs',[]);
P = impose_default_value(P,'maxreads',10000);
P = impose_default_value(P,'Wr',250);                % total window radius
P = impose_default_value(P,'Nr',25);                 % N window radius
P = impose_default_value(P,'mask_central_region',150);
P = impose_default_value(P,'cheat',false);

targleft = round((2*P.Wr - P.mask_central_region)/2);
if targleft<1, error('2*P.Wr < P.mask_central_region'); end
mask = [num2str(targleft) ',' num2str(P.mask_central_region)];

if ~isempty(P.SNP_mask_BAMs), require_fields(P.SNP_mask_BAMs,{'individual','bam'}); end

flds1 = {'individual','name','leftprimer','rightprimer'};
flds2 = {'chr1','pos1','str1','chr2','pos2','str2'};
require_fields(X,union(flds1,flds2));
X = make_numeric(X,flds2);
nx = slength(X);
if length(unique(X.name))<nx, error('X.name must be unique'); end

% build validation target sequences

fprintf('Building target sequences... '); tic

if isempty(P.SNP_mask_BAMs), step = 100; else step = 1; end

X.target = cell(nx,1);
for i=1:nx, if ~mod(i,step), fprintf('%d/%d ',i,nx); end
  dA = masked_genome_region(X.individual{i},X.chr1(i),X.pos1(i)-P.Wr,X.pos1(i)+P.Wr-1);
  if X.str1(i)==1, dA = rc(dA); end
  dB = masked_genome_region(X.individual{i},X.chr2(i),X.pos2(i)-P.Wr,X.pos2(i)+P.Wr-1);
  if X.str2(i)==0, dB = rc(dB); end
  X.target{i} = [dA(1:P.Wr-P.Nr) repmat('N',1,P.Nr*2) dB(P.Wr+P.Nr+1:end)];
end, fprintf('\n');

toc

  function dna = masked_genome_region(individual,chr,st,en)
    if ~exist('en','var'), disp(88); keyboard; end
    dna = genome_region(chr,st,en);
    bases = 'ACGT';
    if ~isempty(P.SNP_mask_BAMs)
      idx = find(strcmp(individual,P.SNP_mask_BAMs.individual));
      B = cell(length(idx),1);
      for bi=1:length(idx)
        [tmp B{bi}] = pull_from_bam(P.SNP_mask_BAMs.bam{idx(bi)},chr,st,en,P);
      end
      B = cat(1,B{~cellfun('isempty',B)});
      range = (st:en)';
      B = B(B(:,2)>25,:);   % keep only high-quality bases
      ntot = as_row(histc(B(:,4),range));
      B = B(B(:,1)>=64,:);  % keep only nonref ACGT
      nmut = as_row(histc(B(:,4),range));
      fracmut = nmut./ntot;
      % number of unique nonref bases at each position
      nnrb = zeros(1,length(range));
      uniquenonrefbase = nan(1,length(range));
      hasdel = false(1,length(range));
      for k=range(1):range(end),j=k-range(1)+1;
        q = unique(B(B(:,4)==k,1));
        nnrb(j) = length(q);
        if length(q)==1, uniquenonrefbase(j) = q; end
        hasdel(j) = ismember(-100,q);
      end
      hznref = (fracmut>=0.9 & ~hasdel & nnrb==1);   % homozygous SNP in tumor+normal
      ismut = find(ntot>=10 & fracmut>=0.15);
      fprintf('\n');
      disp([ismut' ntot(ismut)' nmut(ismut)' nnrb(ismut)' hasdel(ismut)' uniquenonrefbase(ismut)' hznref(ismut)'])
      for k=1:length(ismut),j=ismut(k);
        if hznref(j)
          nonrefbase = bases(uniquenonrefbase(j)-64);
        else
          nonrefbase = 'N';
        end
        dna(j) = nonrefbase;
      end
    end
  end

% ALIGN PRIMERS TO TARGET SEQUENCE

X.check_amplicon = repmat({''},nx,1);
X.check_result = repmat({''},nx,1);


scmat = -ones(size(nuc44)); for z=1:4, scmat(z,z)=1; end
swparams = {'gapopen',20,'alphabet','nt','scoringmatrix',scmat};

for i=1:nx
  msg = '';
  target = X.target{i};
  L.primerseq = X.leftprimer{i};
  L.primerlen = length(L.primerseq);
  try
    [L.a L.b L.c] = swalign(L.primerseq,target,swparams{:});
    L.matchpos = L.c(2);
    L.matchlen = size(L.b,2);
  catch me
    msg = [msg '  LEFT_ALIGN_EXCEPTION'];
    L.matchpos = nan;
    L.matchlen = 0;
  end
  if (P.cheat), target(1:P.Wr) = repmat('N',1,P.Wr); end
  R.primerseq = X.rightprimer{i};
  R.primerlen = length(R.primerseq);
  try
    [R.a R.b R.c] = swalign(rc(R.primerseq),target,swparams{:});
    R.matchpos = R.c(2);
    R.matchlen = size(R.b,2);
  catch me
    msg = [msg '  RIGHT_ALIGN_EXCEPTION'];
    R.matchpos = nan;
    R.matchlen = 0;
  end
  if (L.primerlen-L.matchlen>10), msg = [msg '  LEFT_PRIMER_BAD']; L.matchpos=nan; end
  if (R.primerlen-R.matchlen>10), msg = [msg '  RIGHT_PRIMER_BAD']; R.matchpos=nan; end
  X.ampstart(i,1) = L.matchpos;
  X.ampend(i,1) = R.matchpos+R.matchlen-1;
  X.amplen(i,1) = X.ampend(i)-X.ampstart(i)+1;
  X.amplicon{i,1} = X.target{i,1}(X.ampstart(i):X.ampend(i));
  X.check_amplicon{i,1} = sprintf('amplicon %d-%d (%d)',X.ampstart(i),X.ampend(i),X.amplen(i));
  if (X.amplen(i)<200), msg = [msg '  AMPLICON_TOO_SHORT']; end
  if (X.amplen(i)>400), msg = [msg '  AMPLICON_TOO_LONG']; end
  if isempty(msg), msg = 'OK'; end
  X.check_result{i,1} = msg;
end

fprintf('\nAmplicon lengths:  min %d   max %d  median %d  std %d\n',...
  nanmin(X.amplen),nanmax(X.amplen),nanmedian(X.amplen),nanstd(X.amplen));

fprintf('\nCheck results:');
count(X.check_result);










end % main function
