function C = load_bygene_total_and_maf(individual_id, indir,varargin)
if nargin < 3
  categdir = '/xchip/cga1/annotation/db/tracks/hg19/gc29zc65';
  all_newbases = false; 
else 
  categdir = varargin{1};
  all_newbases = varargin{2};
end 

if isempty(indir) 
  bygene = ['/xchip/cga1/lawrence/mut/20101224/wgs/wig/' individual_id '.coverage.bygene.mat'];
  total = ['/xchip/cga1/lawrence/mut/20101224/wgs/wig/' individual_id '.coverage.total.mat'];
  maf = ['/xchip/cga1/lawrence/mut/20101224/wgs/maf/' individual_id '.gcsz29p.maf.mat'];
else 
  bygene = fullfile(indir, [individual_id '.coverage.bygene.mat']);
  total = fullfile(indir , [individual_id '.coverage.total.mat']);
  maf = fullfile(indir, [individual_id '.maf.annotated']);
  maf_tumor = fullfile(indir, [individual_id '-Tumor.maf.annotated']);
end

try 
  demand_files({bygene;total;maf});
catch 
  demand_files({bygene;total;maf_tumor});
  maf = maf_tumor
end 

% load coverage
fprintf('Loading coverage...\n');
load(bygene,'C')
C = rename_field(C,'cov','N');
if contains(categdir, 'gc29zc65nb') 
keyboard
end 
%keyboard
C.cat = parse_in(C.cat,'name','(good|bad):(nonconserved|conserved):(IGR|intron|exon|UTR):([ACGT]) in ([ACGT])_([ACGT])',...
   {'goodbad','cons','zone','base','left','right'});
load(total,'N')
C.N.c = N; clear N;

% load mutations
fprintf('Loading mutations...\n');
if isempty(indir)
  load(maf,'M')
  M = rename_field(M,'gcsz29p','gc29zc65');
else 
  M = load_struct(maf); 
%  keyboard
  M = move_to_simple_fieldnames(M);
  M = make_numeric(M, {'start'});
  M = reorder_struct(M, ~strcmp(M.chr, 'M') & ~strcmp(M.chr, 'chrM'));
  M.chr = convert_chr(M.chr);
  M.gc29zc65 = get_context(M.chr, M.start, categdir);
end

%keyboard
%M = rename_field(M,'gcsz29p','gc29zc65');   % (forgot to rename the field when we saved it)
M.gidx = listmap(M.gene,C.gene.name);
%M.gidx = mapacross(M.gene,C.gene.name,1:slength(C.gene));   % (equivalent)
%keyboard
C.mut = M; %clear M;
if all_newbases 
  C.n.gc = nan(slength(C.gene), slength(C.cat), 5);
  ma = reorder_struct(M, strcmp(M.tum_allele2, 'A'));
  mc = reorder_struct(M, strcmp(M.tum_allele2, 'C'));
  mg = reorder_struct(M, strcmp(M.tum_allele2, 'G'));
  mt = reorder_struct(M, strcmp(M.tum_allele2, 'T'));
  mn = reorder_struct(M, ~strcmp(M.tum_allele2, 'A') & ~strcmp(M.tum_allele2, 'C') & ~strcmp(M.tum_allele2, 'G') & ~strcmp(M.tum_allele2, 'T'));
  C.n.gc(:, :, 1) = hist2d_fast(ma.gidx,ma.gc29zc65,1,slength(C.gene),1,slength(C.cat));
  C.n.gc(:, :, 2) = hist2d_fast(mc.gidx,mc.gc29zc65,1,slength(C.gene),1,slength(C.cat));
  C.n.gc(:, :, 3) = hist2d_fast(mg.gidx,mg.gc29zc65,1,slength(C.gene),1,slength(C.cat));
  C.n.gc(:, :, 4) = hist2d_fast(mt.gidx,mt.gc29zc65,1,slength(C.gene),1,slength(C.cat));
  C.n.gc(:, :, 5) = hist2d_fast(mn.gidx,mn.gc29zc65,1,slength(C.gene),1,slength(C.cat));
%  keyboard
else 
  C.n.gc = hist2d_fast(C.mut.gidx,C.mut.gc29zc65,1,slength(C.gene),1,slength(C.cat));
end 
C.n.c = histc(C.mut.gc29zc65,1:slength(C.cat));
