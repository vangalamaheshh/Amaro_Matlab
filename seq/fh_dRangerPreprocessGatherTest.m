function fh_dRangerPreprocessGatherTest(outdir)
% fh_dRangerPreprocessGather(outdir)
%
% outdir = <SAMPLE_ID>.dRanger_input, will be name of output directory
%
% (1) Creates output directory
% (2) Moves output of all scatter jobs to the output directory
% (3) Loads, joins, and de-dups chr*.weird, saves as dRanger_input.mat
% (4) Keeps all filtering-metric files (*.nuwp and *.fmapqz) as is
% (5) Collapses *.isz files by chromsome and agglomerates by library
% (6) Deletes intermediate results
%
% Mike Lawrence 2010-01-22

fprintf('fh_dRangerPreprocessGather\n');
fprintf(['  outdir = "' outdir '"\n']);

% move results of all scatter jobs to the output directory

if ~exist(outdir,'dir'), mkdir(outdir); end
status = system(['mv scatter*/chr* ' outdir]);
if status~=0, fprintf('Warning: file move failed\n'); end

cd(outdir);

% load and join weird reads

dirchr=dir(['chr*.weird']);
NC=length(dirchr);
Xj=cell(NC,1);
Xu=cell(NC,1);
for c=1:NC
  fprintf('CHR%d\n',c);
  fprintf('   Loading... ');tic
  fname = ['chr' num2str(c) '.weird'];
  X = read_table(fname,repmat('%f',1,13),char(9),0);toc
  X = cat(2,X.dat{:});
  fprintf('   Joining intrachromosomal... ');tic
  [Xj{c} Xu{c}] = join_matching_rows(X);toc
end

fprintf('Concatenating result... ');tic
Xj = cat(1,Xj{:});
Xu = cat(1,Xu{:});
toc

fprintf('Joining interchromosomal... ');tic
[Xa Xb] = join_matching_rows(Xu);
clear Xu;
X = cat(1,Xj,Xa,Xb);
clear Xa Xb;
toc

% starting order:   rgrp namenumber chr1 start1 end1 strand1 qual1 chr2 start2 end2 strand2 qual2 flip
% desired order: namenumber chr1 strand1 start1 end1 chr2 strand2 start2 end2 qual1 qual2 rgrp
X = X(:,[2 3 6 4 5 8 11 9 10 7 12 1]);

% de-dup
tot = size(X,1);
tic,fprintf('De-duping... ');
[u ui uj] = unique(X(:,2:9),'rows');
X = X(ui,:);
dd = length(ui);
toc,fprintf('Number of records:\n\tTotal: %d\n\tDe-duped: %d\n',tot,dd);

% save
fprintf('Saving... ');tic
save('all.weird.pairs.mat','X','-v7.3');toc

if exist('all.weird.pairs.mat','file')
  system('rm chr*.weird');
else
  error('error saving all.weird.pairs.mat');
end

% collapse isz files by chromosome

A = []; I = [];
for c = 1:NC
  tmp = load_matrix(['chr' num2str(c) '.isz']);
  if c==1, A = tmp(:,1); I = tmp(:,2:end); else I(:,:,c) = tmp(:,2:end); end
end
I = sum(I,3);  % collapse chromosomes

out = fopen('all.isz', 'wt');
for i=1:size(I,1);
  fprintf(out,'%d', A(i));
  for j=1:size(I,2);
    fprintf(out, '\t%d', I(i,j));
  end
  fprintf(out,'\n');
end
fclose(out);

if exist('all.isz','file')
  system('rm chr*.isz');
else
  error('Failed in writing all.isz');
end

% agglomerate by library (for dRangerLocal)
% (keep unagglomerated version)

agglom_iszfile('all.isz','all.agglom.isz');

% create links to the most important files

cd('..');
sample_id = regexprep(outdir,'\.dRanger_input','');
force_link([outdir '/all.isz'],[sample_id '.all.isz']);
force_link([outdir '/all.agglom.isz'],[sample_id '.all.agglom.isz']);
force_link([outdir '/all.weird.pairs.mat'],[sample_id '.all.weird.pairs.mat']);

% done!
