function fh_InsertSizeByLaneGather(outputfile,flag)
% fh_InsertSizeByLaneGather(outputfile)
%
% produces a single output file which is a combination of all the chr*.isz files
%
% Mike Lawrence / Doug Voet 2009

if ~exist('outputfile','var'), error('input parameter required: outputfile location'); end
if ~exist('flag','var'), flag = false; end

% GATHER
A = []; I = [];
for c = 1:24
  if flag  % for testing purposes
    tmp = load_matrix(['chr' num2str(c) '.isz']);
  else     % for normal Firehose use
    tmp = load_matrix([sprintf('scatter.%010d/chr%d.isz', c, c)]);
  end
  if c==1, A = tmp(:,1); I = tmp(:,2:end); else I(:,:,c) = tmp(:,2:end); end
end
I = sum(I,3);  % collapse chromosomes

out = fopen(outputfile, 'wt');
for i=1:size(I,1);
  fprintf(out,'%d', A(i));
  for j=1:size(I,2);
    fprintf(out, '\t%d', I(i,j));
  end
  fprintf(out,'\n');
end
fclose(out);

if ~flag
  for c=1:24, delete([sprintf('scatter.%010d/chr%d.isz', c, c)]); end
end

fprintf('Finished successfully.\n');
