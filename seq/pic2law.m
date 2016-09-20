function pic2law(listfile)

T = load_struct(listfile);
require_fields(T,{'squid','tumcode','normcode','sample','tumver','normver'});
%e.g.   squid =  G2095
%     tumcode = V0D56W
%    normcode = V0D56X
%      sample = mm/0319/wgs
%      tumver = v3
%     normver = v2

bd1 = '/seq/picard/by_project';
bd2 = '/xchip/tcga_scratch/lawrence';

for i=1:slength(T)
  fprintf('\nPROCESSING SAMPLE %s\n',T.sample{i});
  dir1 = [bd1 '/' T.squid{i}];
  dir2 = [bd2 '/' T.sample{i}];
  if ~exist(dir1,'dir'), error('Not found: %s',dir1); end
  if ~exist(dir2,'dir')
    fprintf('Creating: %s\n',dir2);
    mkdir(dir2);
  end
  tdir = [dir1 '/' T.tumcode{i} '/' T.tumver{i}];
  force_link([tdir '/' T.tumcode{i} '.bam'], [dir2 '/tumor.bam'] );
  force_link([tdir '/' T.tumcode{i} '.bai'], [dir2 '/tumor.bai'] );
  force_link([dir2 '/tumor.bam'],[dir2 '/tumor.bam.bai']);
  ndir = [dir1 '/' T.normcode{i} '/' T.normver{i}];
  force_link([ndir '/' T.normcode{i} '.bam'], [dir2 '/normal.bam'] );
  force_link([ndir '/' T.normcode{i} '.bai'], [dir2 '/normal.bai'] );
  force_link([dir2 '/normal.bam'],[dir2 '/normal.bam.bai']);
end

