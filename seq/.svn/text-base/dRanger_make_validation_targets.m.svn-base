function dRanger_make_validation_targets(X,P)
% Mike Lawrence 2009-10-26

if ~exist('P','var'), P = []; end

if ~isfield(P,'build')
  fprintf('Assuming hg18\n');
  P = impose_default_value(P,'build','hg18');
end
P = impose_default_value(P,'Wr',800);                % total window radius
P = impose_default_value(P,'Nr',25);                 % N window radius
P = impose_default_value(P,'firstno',1);             % numbering scheme start
P = impose_default_value(P,'outstem','*required*');  % numbering scheme start
P = impose_default_value(P,'namestem','R');
P = impose_default_value(P,'SNP_mask_BAMs',[]);

flds1 = {'individual','name'};
flds2 = {'chr1','pos1','str1','chr2','pos2','str2'};
require_fields(X,union(flds1,flds2));
X = make_numeric(X,flds2);
nx = slength(X);
if length(unique(X.name))<nx, error('X.name must be unique'); end

FA=[]; FB=[]; FT=[];
for i=1:slength(X), fprintf('%d/%d\n',i,slength(X));
  dA = masked_genome_region(X.individual{i},X.chr1(i),X.pos1(i)-P.Wr,X.pos1(i)+P.Wr-1,P.build);
  if X.str1(i)==1, dA = rc(dA); end
  dB = masked_genome_region(X.individual{i},X.chr2(i),X.pos2(i)-P.Wr,X.pos2(i)+P.Wr-1,P.build);
  if X.str2(i)==0, dB = rc(dB); end
  dT = [dA(1:P.Wr-P.Nr) repmat('N',1,P.Nr*2) dB(P.Wr+P.Nr+1:end)];
%  seqname{i,1} = [P.namestem num2str(i+P.firstno-1)];
  seqname{i,1} = X.name{i};
  FA.header{i,1} = [seqname{i} '_A']; FA.seq{i,1} = dA;
  FB.header{i,1} = [seqname{i} '_B']; FB.seq{i,1} = dB;
  FT.header{i,1} = [seqname{i} '_T']; FT.seq{i,1} = dT;
%find(dA=='N'),find(dB=='N')
end

name = P.outstem;
save_fasta(FA,[name '_A.fasta']);
save_fasta(FB,[name '_B.fasta']);
save_fasta(FT,[name '_T.fasta']);

K = [];
K.seqname = seqname;
K = merge_structs({K,X});
save_struct(K,[name '_rearrs.txt']);

  function dna = masked_genome_region(individual,chr,st,en,build)
    dna = genome_region(chr,st,en,build);
    bases = 'ACGT';
    if ~isempty(P.SNP_mask_BAMs)
      idx = find(strcmp(individual,P.SNP_mask_BAMs.individual));
      B = cell(length(idx),1);
      if isempty(idx), fprintf('No BAMS match %s\n', individual); return; end
      for bi=1:length(idx)
        [tmp B{bi}] = pull_from_bam(P.SNP_mask_BAMs.bam{idx(bi)},chr,st,en,P);
      end
      B = cat(1,B{~cellfun('isempty',B)});
      if isempty(B), fprintf('BAMS contain no data for this region'); return; end
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

%  function dna = masked_genome_region(chr,st,en)
%    dna = genome_region(chr,st,en);
%    if ~isempty(P.SNP_mask_BAM)
%      for b=1:length(P.SNP_mask_BAM)
%        [R B S] = pull_from_bam(P.SNP_mask_BAM{b},chr,st,en);
%        range = (st:en)';
%        B = B(B(:,2)>25,:);   % keep only high-quality bases
%        ntot = as_row(histc(B(:,4),range));
%        B = B(B(:,1)>=64,:);  % keep only nonref ACGT
%        nmut = as_row(histc(B(:,4),range));
%        fracmut = nmut./ntot;
%        ismut = find(ntot>=10 & fracmut>=0.2);
%        dna(ismut) = 'N';
%      end
%    end
%  end

end % main function
