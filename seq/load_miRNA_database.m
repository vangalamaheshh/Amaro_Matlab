function R = load_miRNA_database(build)
%
% loads miRNA databases from stored file
%

if ~exist('build','var'), build = 'hg18'; end
dirname = ['/xchip/tcga/gbm/analysis/lawrence/genome/' build '/miRNA/'];
fname = [dirname 'hsa_gff.txt'];

fprintf('Loading miRNA database...\n');

% ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/hsa.gff
% remove all header lines and save as hsa_gff.txt

X = load_struct(fname,'%f%*s%s%f%f%*s%s%*s%s',0);

R = [];
R.type = X.col2;
R.chr = X.col1;
R.start = X.col3;
R.end = X.col4;
R.strand = X.col5;
tmp_acc =  regexp(X.col6,'ACC="([^"]*)"','tokens');
tmp_id =  regexp(X.col6,'ID="([^"]*)"','tokens');
R.acc = cell(slength(R),1);
R.id = cell(slength(R),1);
for i=1:slength(X)
  R.acc(i) = tmp_acc{i}{1};
  R.id(i) = tmp_id{i}{1};
end
