function X = load_perlane_chunk_data(samples)
% Mike Lawrence 2009-07-17
% modified for more general use 2009-08-12

basedir = '/xchip/tcga_scratch/lawrence';
X = []; X.dat = []; X.lane = []; X.chunk = []; lane = 1;
fprintf('Loading samples: ');
for i=1:length(samples), fprintf('%d ',i);
  m = sortrows(load_matrix([basedir '/' samples{i} '/chunk10Mb_cov_per_lane_normal.txt']));
  if i==1, X.chunk.num = m(:,1);X.chunk.chr = m(:,2);X.chunk.start = m(:,3);X.chunk.end = m(:,4);end
  lt = []; ltidx = 1;
  for j=1:size(m,2)-4
    X.dat(:,lane) = m(:,4+j);
    lt.sample{ltidx,1} = samples{i};
    lt.tn{ltidx,1} = 'N';
    lt.short{ltidx,1} = regexprep(samples{i},'^[^/]+/([^/]+)/[^/]+$','$1');
    lt.id(ltidx,1) = j-1;
    ltidx=ltidx+1;
    lane=lane+1;
  end
  fname = [basedir '/' samples{i} '/normal.bam.lanetable'];
  if exist(fname,'file')
    lt = merge_structs({lt,load_lanetable(fname)});
  end
  if i==1, X.lane = lt; else X.lane = concat_structs({X.lane,lt}); end
  m = sortrows(load_matrix([basedir '/' samples{i} '/chunk10Mb_cov_per_lane_tumor.txt']));
  lt = []; ltidx = 1;
  for j=1:size(m,2)-4
    X.dat(:,lane) = m(:,4+j);
    lt.sample{ltidx,1} = samples{i};
    lt.tn{ltidx,1} = 'T';
    lt.short{ltidx,1} = regexprep(samples{i},'^[^/]+/([^/]+)/[^/]+$','$1');
    lt.id(ltidx,1) = j-1;
    ltidx=ltidx+1;
    lane=lane+1;
  end
  fname = [basedir '/' samples{i} '/tumor.bam.lanetable'];
  if exist(fname,'file')
    lt = merge_structs({lt,load_lanetable(fname)});
  end
  X.lane = concat_structs({X.lane,lt});
end, fprintf('\n');
[X.nchunks, X.nlanes] = size(X.dat);


