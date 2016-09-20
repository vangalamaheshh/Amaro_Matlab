function C = sum_binary_coverage_and_collapse_by_gene(categdir,targlist,binfiles)

C = [];
C.file = [];
C.file.categdir = categdir;
C.file.targlist = targlist;
C.file.binfiles = binfiles;

% get number of categories
C.cat = load_struct([categdir '/categs.txt']);
if strcmp(C.cat.num{1},'0'), error('Category list includes a zero category!'); end
C.ncat = slength(C.cat);

% load target file
C.targ = load_target_file(targlist);
demand_fields(C.targ,{'gene'});
C.nt = slength(C.targ);

% collapse targets to genes
[C.gene.name tmp C.targ.gidx] = unique(C.targ.gene);
C.ng = slength(C.gene);

% allocate storage for result
C.cov = zeros(C.ng,C.ncat,'int32');

% check input files
nf = length(binfiles);
demand_files(binfiles);
ct = C.nt*C.ncat;
fstring = cell(nf,1);
for f=1:nf
  d = dir(binfiles{f});
  sz = d.bytes;
  if sz==ct, fstring{f}='int8=>int32';
  elseif sz==ct*2, fstring{f}='int16=>int32';
  elseif sz==ct*4, fstring{f}='int32=>int32';
  else error('size mismatch: %s',fname);
  end
end
endianness = 'b';

% load data
for f=1:nf
  fprintf('%d/%d %s\n',f,nf,binfiles{f});
  fh = fopen(binfiles{f});
  for i=1:C.nt, if ~mod(i,1000), fprintf('%d/%d ',i,C.nt); end
    C.cov(C.targ.gidx(i),:) = C.cov(C.targ.gidx(i),:) + fread(fh,C.ncat,fstring{f},endianness)';
  end, fprintf('\n');
  fclose(fh);
end
