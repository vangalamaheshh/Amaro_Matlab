function C = load_mutation_spectrum_raw_data(indiv)

dr = '/xchip/cga1/lawrence/mut/20101224/wgs';
fprintf('Loading per-gene coverage\n');
load([dr '/wig/' indiv '.coverage.bygene.mat'],'C');
fprintf('Loading total coverage\n');
load([dr '/wig/' indiv '.coverage.total.mat'],'N');
fprintf('Loading mutations\n');
load([dr '/maf/' indiv '.gcsz29p.maf.mat'],'M');
C.gene.name{end+1} = '__IGR__';
C.ng = C.ng + 1;
C.cov.gc(end+1,:) = N'-sum(C.cov.gc);
M.gidx(isnan(M.gidx)) = C.ng; % IGR mutations
C.mut = M;
clear M N
C.n.gcn = nan(C.ng,C.ncat,4);
for i=1:4
  idx = find(C.mut.newbaseidx==i);
  C.n.gcn(:,:,i) = hist2d_fast(C.mut.gidx(idx),C.mut.gcsz29p(idx),1,C.ng,1,C.ncat);
end

% add strand and footprint to gene
G=get_refseq_introns_and_exons('hg18');
idx = listmap(C.gene.name,G.name);
C.gene.chr = nansub(G.chr,idx);
C.gene.strand = nansub(G.strand,idx);
C.gene.footprint = nansub(G.footprint,idx);

% add expression (using CRC normals--should download a MEL expr set from GEO)
L = load_struct('/xchip/cga1/lawrence/mut/20110103/expr/crc_GDS2947.txt','%s%f%f%f');
C.gene.expr = mapacross(C.gene.name,L.name,L.medlog_tumor);



