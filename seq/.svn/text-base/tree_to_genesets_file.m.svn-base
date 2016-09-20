function tree_to_genesets_file(genenames,tree,outfile,num)
% tree_to_genesets_file(genenames,tree,outfile,num)
%
% Mike Lawrence 2010-03-08

nc = size(tree,1);
ng = length(genenames);

if ~exist('num','var'), num = nc; end

fprintf('Computing subsets: ');
X = [eye(ng);zeros(nc,ng)];
for i=1:nc, if ~mod(i,1000), fprintf('%d/%d ',i,nc); end
  X(ng+i,:) = X(tree(i,1),:) + X(tree(i,2),:);
end,fprintf('\n');

fprintf('Writing file: ');
f = fopen(outfile,'wt');
for i=1:num, if ~mod(i,1000), fprintf('%d/%d ',i,nc); end
  fprintf(f,'Cluster%d\tdist=%0.3f',i,tree(i,3));
  fprintf(f,'\t%s',genenames{X(ng+i,:)>0});
  fprintf(f,'\n');
end,fprintf('\n');
fclose(f);



