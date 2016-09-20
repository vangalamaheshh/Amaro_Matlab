% download chr_rpt files using "download" script


% 2009-10-07
% human build "9606" is dbSNP release 130

%set baseurl = 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/chr_rpts';
%set basedir = '/xchip/tcga_scratch/lawrence/db/dbsnp/130';

%cd $basedir;
%foreach chr (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
%  bsub wget $baseurl/chr_$chr.txt.gz
%end

%cd $basedir
%foreach chr (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
%  gunzip chr_$chr.txt.gz
%end

% manually constructed file "header" based on lines 4+5 of the textfiles

%cd $basedir
%foreach chr (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
%  tail +8 chr_$chr.txt | grep 'reference$' > body
%  cat header body > chr_$chr.edited.txt
%end


% convert textfiles to binary "*.dbSNP" files
% which are arrays of integers
%   each nonzero element tells the rs of the dbSNP reported at that position

basedir = '/xchip/tcga_scratch/lawrence/db/dbsnp/130';

java_classpath = [...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

len = load_chrlen;
for i=1:24
  if i<23, c = num2str(i);elseif i==23, c='X';else c='Y'; end
  infile = [basedir '/chr_' c '.txt'];
  outfile = [basedir '/chr' num2str(i) '.dbSNP'];
  cmd = ['"java -Xmx2g -classpath ' java_classpath ' '...
           'Convert_dbSNP ' infile ' ' outfile ' ' c ' ' num2str(len(i)) '"'];
  banner = ['dbSNP_' c];
  bsub(cmd,banner);
end
