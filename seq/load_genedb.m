function D = load_genedb(build)
% load_genedb(build)
% 
% loads combined gene database from stored files
%
% Currently includes: RefSeq, CCDS, Ensembl
%
% adjusts all starting coordinates forward by one nucleotide 
% to undo confusing refseq convention of zero-based start, one-based end
%

if ~exist('build','var'), build = 'hg18'; end
dirname = ['/xchip/tcga/gbm/analysis/lawrence/genome/' build '/'];
matname = [dirname 'D.mat'];

fprintf('Loading gene database...\n');
if strcmp(build,'hg18')

  if exist(matname,'file')
    load(matname,'D');

  else
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/ccdsGene.txt.gz
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/ensGene.txt.gz

    D1 = load_struct([dirname 'ccdsGene.txt']);
    D2 = load_struct([dirname 'ensGene.txt']);
    D3 = load_struct([dirname 'refGene.txt']);
    D1 = rmfield(D1,'score');
    D2 = rmfield(D2,'score');
    D3 = rmfield(D3,'id');
    D1.db = repmat({'CCDS'},slength(D1),1);
    D2.db = repmat({'Ensembl'},slength(D2),1);
    D3.db = repmat({'RefSeq'},slength(D3),1);
    D = combine_structs({D1,D2,D3});
    fprintf('Converting...\n');
    D.txStart = str2double(D.txStart) + 1;
    D.txEnd = str2double(D.txEnd);
    D.cdsStart = str2double(D.cdsStart) + 1;
    D.cdsEnd = str2double(D.cdsEnd);
    D.exonCount = str2double(D.exonCount);
    for i=1:slength(D)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(D)); end
      D.exonStarts{i} = str2double(split(D.exonStarts{i}(1:end-1),','))+1;
      D.exonEnds{i} = str2double(split(D.exonEnds{i}(1:end-1),','));
      D.exonFrames{i} = str2double(split(D.exonFrames{i}(1:end-1),','));
    end

    save(matname,'D');
  end

elseif strcmp(build,'hg17')
  if exist(matname,'file')
    load(matname,'D');

  else
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/refGene.txt.gz
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/ccdsGene.txt.gz
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/ensGene.txt.gz

    T = read_table([dirname 'refGene.txt'],'%s%s%s%f%f%f%f%f%s%s',...
       char(9),0,'whitespace','\b\r');
    [D1.name,D1.chrom,D1.strand,D1.txStart,D1.txEnd,D1.cdsStart,D1.cdsEnd,...
      D1.exonCount,D1.exonStarts,D1.exonEnds]  = deal(T.dat{:});
    T = read_table([dirname 'ensGene.txt'],'%s%s%s%f%f%f%f%f%s%s',...
       char(9),0,'whitespace','\b\r');
    [D2.name,D2.chrom,D2.strand,D2.txStart,D2.txEnd,D2.cdsStart,D2.cdsEnd,...
      D2.exonCount,D2.exonStarts,D2.exonEnds]  = deal(T.dat{:});
    T = read_table([dirname 'ccdsGene.txt'],'%*f%s%s%s%f%f%f%f%f%s%s%*s%*s%*s%*s%*s',...
       char(9),0,'whitespace','\b\r');
    [D3.name,D3.chrom,D3.strand,D3.txStart,D3.txEnd,D3.cdsStart,D3.cdsEnd,...
      D3.exonCount,D3.exonStarts,D3.exonEnds]  = deal(T.dat{:});

    D3.db = repmat({'CCDS'},slength(D3),1);
    D2.db = repmat({'Ensembl'},slength(D2),1);
    D1.db = repmat({'RefSeq'},slength(D1),1);
    D = combine_structs({D1,D2,D3});
    fprintf('Converting...\n');
    D.txStart = D.txStart + 1;
    D.cdsStart = D.cdsStart + 1;
    for i=1:slength(D)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(D)); end
      D.exonStarts{i} = str2double(split(D.exonStarts{i}(1:end-1),','))+1;
      D.exonEnds{i} = str2double(split(D.exonEnds{i}(1:end-1),','));
    end

    save(matname,'D');
  end

else error('Unsupported build %s',build);
end

% filtering of dubious Ensembl transcripts
%    ENST00000380949    8.9 Mb
%    ENST00000402964    4.3 Mb
%    ENST00000406903    50.9 Mb
%    ENST00000402654    4.0 Mb
%    ENST00000402917    10.6 Mb

D.txlen = D.txEnd - D.txStart + 1;
D = reorder_struct(D,D.txlen<3e6);

