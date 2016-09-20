infile = '/xchip/cga1/lawrence/capture/whole_exome_agilent_designed_120.targets.interval_list';
outfile = ['/xchip/cga1/lawrence/capture/' ...
           'whole_exome_agilent_designed_120.targets.interval_list.GENESGC_noheader.txt'];
process_targets_interval_list(infile,outfile);

fc2k = '/xchip/cga1/lawrence/capture/cancer_2000gene_shift170.targets.interval_list_NUMSONLY.txt';
fc6k = '/xchip/cga1/lawrence/capture/tcga_6k_genes.targets.interval_list_NUMSONLY.txt';
fwe = ['/xchip/cga1/lawrence/capture/' ...
           'whole_exome_agilent_designed_120.targets.interval_list.GENESGC_noheader.txt'];

c2k = rename_fields(load_struct(fc2k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = rename_fields(load_struct(fc6k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = reorder_struct(c6k,~isnan(c6k.chr));
we = rename_fields(load_struct(fwe,'%s%f%f%f%f',0),colx(1:5),{'gene','chr','start','end','gc'});
m = understand_nested_target_sets({we,c6k,c2k});
save('/xchip/cga1/lawrence/capture/c2k_c6k_we16k_nesting.mat','m');

fc2k = '/xchip/cga1/lawrence/capture/cancer_2000gene_shift170.targets.interval_list_NUMSONLY.txt';
fc6k = '/xchip/cga1/lawrence/capture/tcga_6k_genes.targets.interval_list_NUMSONLY.txt';
c2k = rename_fields(load_struct(fc2k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = rename_fields(load_struct(fc6k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = reorder_struct(c6k,c6k.chr>=1 & c6k.chr<=24);
m = understand_nested_target_sets({c6k,c2k});
save('/xchip/cga1/lawrence/capture/c2k_c6k_nesting.mat','m');



%%  make wiggle file of c6k targets
%% 
fc6k = '/xchip/cga1/lawrence/capture/tcga_6k_genes.targets.interval_list_NUMSONLY.txt';
c6k = rename_fields(load_struct(fc6k,'%f%f%f',0),colx(1:3),{'chr','start','end'});
c6k = reorder_struct(c6k,~isnan(c6k.chr));
g = Genome();
for i=1:slength(c6k)%, if ~mod(i,100), fprintf('%d/%d ',i,slength(c6k)); end
  g.setContents(c6k.chr(i),c6k.start(i),c6k.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(c6k),g.totCount());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/tcga_6k_genes.targets.wig');

%%  make wiggle file of weRefseq targets
%%
f = '/xchip/cga1/lawrence/capture/whole_exome_refseq_coding.targets.interval_list.GENESGC.txt';
t = rename_fields(load_struct(f,'%s%f%f%f%f',0),colx(1:5),{'gene','chr','start','end','gc'});
t = reorder_struct(t,~isnan(t.chr));
g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/whole_exome_refseq_coding.targets.wig');


%%  2010-02-17
%%  try matching capture target sets directly to Refseq transcripts

R = load_refseq;
match_targets_to_refseq('/xchip/cga1/lawrence/capture/tcga_6k_genes.targets.interval_list',R);

%%  2010-06-01

outfile = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_dec2009.txt';
prepare_mutsig_target_list(outfile,'hg18',0,0.20);

%%  2010-06-23
%%  make wiggle file of RefSeq exons
%%
f = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_dec2009.txt';
t = rename_fields(load_struct(f,'%s%f%f%f%f%f%f',0),colx(1:7),{'gene','chr','start','end','gc','len','mem'});
t = reorder_struct(t,~isnan(t.chr));
g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_dec2009.wig');

% make reduced version of context65 wiggle
g = Genome();
g.loadWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_dec2009.wig','mask');
g.loadWiggle('/xchip/cga1/lawrence/db/context65/all.wig');
g.saveWiggle('/xchip/cga1/lawrence/db/context65/RefSeq_exons_hg18_dec2009_terr_only.wig');

% make RefSeq exons genelist
G = get_refseq_introns_and_exons('hg18',0);
G = keep_fields(G,{'name'});
save_struct(G,'/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_dec2009.genelist.txt');

%%  2010-08-03
%%  make new RefSeq target list and wiggle files
%%  (1) using newer version of Refseq
%%  (2) including 2-bp splice-site flanks on exons

outfile = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010.txt';
prepare_mutsig_target_list(outfile,'hg18_v2',2,0.20);

f = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010.txt';
t = rename_fields(load_struct(f,'%s%f%f%f%f%f%f',0),colx(1:7),{'gene','chr','start','end','gc','len','mem'});
t = reorder_struct(t,~isnan(t.chr));
clear g; g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010.wig');

% make reduced version of context65 wiggle
clear g; g = Genome();
g.loadWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010.wig','mask');
g.loadWiggle('/xchip/cga1/lawrence/db/context65/all.wig');
g.saveWiggle('/xchip/cga1/lawrence/db/context65/RefSeq_exons_hg18_june2010_terr_only.wig');

% make RefSeq exons genelist
G = get_refseq_introns_and_exons('hg18_v2',2);
% (NOTE: 27 genes have duplicate entries! many on both chrX&chrY)
G = keep_fields(G,{'name'});
save_struct(G,'/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010.genelist.txt');

%%  2010-09-07
%%  fixed bug that was causing 2186 redundant records (overlapping exons due to incomplete collapse)

outfile = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1.txt';
prepare_mutsig_target_list(outfile,'hg18_v2',2,0.20);
% (associated WIG files and genelist can be used as-is)

%%  2010-09-10
%%  recalculate "membership" column on baitset used in MM project

outfile = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_MM.txt';
exome_interval_list = '/seq/references/HybSelOligos/whole_exome_agilent_designed_120/whole_exome_agilent_designed_120.targets.interval_list';
prepare_mutsig_target_list(outfile,'hg18_v2',2,0.20,exome_interval_list);
% (associated WIG files and genelist can be used as-is)

%%  2010-09-14
%%  make hg19 version

build = 'hg19';
outfile = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.txt';
exome_interval_list = ['/xchip/cga/reference/hg19/' ...
                    'whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_100bp_padding_minus_mito' ...
                    '.Homo_sapiens_assembly19.targets.interval_list'];
splicesiteflank = 2;
threshold = 0.20;
prepare_mutsig_target_list(outfile,build,splicesiteflank,threshold,exome_interval_list);

% make mask wiggle
f = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.txt';
t = rename_fields(load_struct(f,'%s%f%f%f%f%f%f',0),colx(1:7),{'gene','chr','start','end','gc','len','mem'});
t = reorder_struct(t,~isnan(t.chr));
clear g; g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.wig');

% make reduced version of context65 wiggle
clear g; g = Genome();
g.loadWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.wig','mask');
g.loadWiggle('/xchip/cga1/lawrence/db/hg19/context65/all.wig');
g.saveWiggle('/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_hg19_june2010_terr_only.wig');

% make RefSeq exons genelist
G = get_refseq_introns_and_exons('hg19',2);
% (NOTE: 26 genes have duplicate entries! many on both chrX&chrY)
G = keep_fields(G,{'name'});
save_struct(G,'/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.genelist.txt');

%%  2010-09-30
%%  for /xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1.txt
%%     add the following 10 genes from WeRefseq
%%           AOF2; FBXL(10,11); JMJD__(1A,1B,1C,2A,2B,2C,3)
%%  save as /xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.txt

% grep 'AOF2|FBXL10|FBXL11|JMJD1A|JMJD1B|JMJD1C|JMJD2A|JMJD2B|JMJD2C|JMJD2D|JMJD3' whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP.txt | sort > plus11.txt
% cat /xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1.txt plus11.txt > /xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.txt

% re-sort by chr-start-end
q = load_struct_noheader('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.txt');
q = make_numeric(q,{'col2','col3','col4'});
q = sort_struct(q,{'col2','col3','col4'});
save_struct_noheader(q,'/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.txt');

% make auxiliary files

% (1) genelist
f = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.txt';
t = rename_fields(load_struct(f,'%s%f%f%f%f%f%f',0),colx(1:7),{'gene','chr','start','end','gc','len','mem'});
t = reorder_struct(t,~isnan(t.chr));
G=[]; G.name = unique(t.gene);
save_struct(G,'/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.genelist.txt');

% (2) wiggle file telling which bases are in the target list
clear g; g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.wig');

% (3) reduced version of context65 wiggle
clear g; g = Genome();
g.loadWiggle('/xchip/cga1/lawrence/capture/RefSeq_exons_hg18_june2010_fix1_plus11.wig','mask');
g.loadWiggle('/xchip/cga1/lawrence/db/context65/all.wig');
g.saveWiggle('/xchip/cga1/lawrence/db/context65/RefSeq_exons_hg18_june2010_fix1_plus11_terr_only.wig');


%%  2010-10-05
%%  make new target lists from Refseq as follows

P = [];
P.refseqbuild = '/xchip/cga1/annotation/db/ucsc/hg18_v2/R.mat';  % (2010-06-19)
P.genomeregionbuild = '/xchip/cga1/annotation/db/ucsc/hg18_v2';
P.conservation_matfile = '/xchip/cga1/lawrence/mm/analysis/20100318_cons/regulatory_regions_29mammal.mat';
P.splicesiteflank = 2;
P.imputed_promoter_length = 3000;
name = 'hg18_june2010_cons29mammal';
outdir = '/xchip/cga1/lawrence/capture'; ensure_dir_exists(outdir);

P.include_introns = true; P.include_UTRs = true; P.include_promoters = true;
P.include_exons = true; P.restrict_to_conserved_regions = false;
save_struct_noheader(generate_target_list(P),[outdir '/' name '.full_transcript.txt']);

P.include_introns = true; P.include_UTRs = true; P.include_promoters = true;
P.include_exons = true; P.restrict_to_conserved_regions = true;
save_struct_noheader(generate_target_list(P),[outdir '/' name '.full_transcript.cons_regions_only.txt']);

P.include_introns = true; P.include_UTRs = true; P.include_promoters = true;
P.include_exons = false; P.restrict_to_conserved_regions = false;
save_struct_noheader(generate_target_list(P),[outdir '/' name '.noncoding.txt']);

P.include_introns = true; P.include_UTRs = true; P.include_promoters = true;
P.include_exons = false; P.restrict_to_conserved_regions = true;
save_struct_noheader(generate_target_list(P),[outdir '/' name '.noncoding.cons_regions_only.txt']);

P.include_introns = false; P.include_UTRs = false; P.include_promoters = false;
P.include_exons = true; P.restrict_to_conserved_regions = false;
save_struct_noheader(generate_target_list(P),[outdir '/' name '.exons.txt']);

P.include_introns = false; P.include_UTRs = false; P.include_promoters = false;
P.include_exons = true; P.restrict_to_conserved_regions = true;
save_struct_noheader(generate_target_list(P),[outdir '/' name '.exons.cons_regions_only.txt']);

%%  2010-12-23
%%  make mask and reduced context65 wigs for hg18 and hg19 versions of common "good exons" list

% hg18
t = load_target_file('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg18.txt');
% mask
clear g; g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg18.mask.wig');
% reduced context65 wig
clear g; g = Genome();
g.loadWiggle('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg18.mask.wig','mask');
g.loadWiggle('/xchip/cga1/lawrence/db/context65/all.wig');
g.saveWiggle('/xchip/cga1/lawrence/db/context65/RefSeq_exons_good_20101221_hg18_terr_only.wig');

% hg19
t = load_target_file('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.txt');
% mask
clear g; g = Genome();
for i=1:slength(t)
  g.setContents(t.chr(i),t.start(i),t.end(i),1);
  fprintf('%d/%d\t%d\n',i,slength(t),g.getTotContents());
end, fprintf('\n');
g.saveWiggle('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.mask.wig');
% reduced context65 wig
clear g; g = Genome();
g.loadWiggle('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.mask.wig','mask');
g.loadWiggle('/xchip/cga1/lawrence/db/hg19/context65/all.wig');
g.saveWiggle('/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_good_20101221_hg19_terr_only.wig');

% hg19 transcripts list
P=[];P.add_promoters = true; % (to match hg18 version)
G=get_refseq_introns_and_exons('hg19',P)
% (promoters weren't actually added!  should have used gene_start+gene_end)
X=keep_fields(G,{'name','chr','tx_start','tx_end'});
X.gc = 0.50*ones(slength(X),1);
X.len = X.tx_end - X.tx_start + 1;
X.mem = ones(slength(X),1);
save_struct_noheader(X,'/xchip/cga1/lawrence/capture/hg19.full_transcript.20110113.txt');


% hg19 transcripts list with 1Kb flanks
P=[];P.add_promoters = false;
G=get_refseq_introns_and_exons('hg19',P)
X=keep_fields(G,{'name','chr','tx_start','tx_end'});
flank = 1000;
X.tx_start = X.tx_start - flank;
X.tx_end = X.tx_end + flank;
X.gc = 0.50*ones(slength(X),1);
X.len = X.tx_end - X.tx_start + 1;
X.mem = ones(slength(X),1);
save_struct_noheader(X,'/xchip/cga1/lawrence/capture/hg19.full_transcript_with_1Kb_flanks.20111020.txt');
% check for overlaps
X = sort_struct(X,{'chr','tx_start'});
X.gap = [nan;X.tx_end(2:end) - X.tx_start(1:end-1)];
sum(X.gap<0) % 23  (the interchromosomal jumps):   perfect!

% hg19 transcripts list with 3Kb and 0kb 5'- and 3'-flanks
% (to match Oncotator)
P=[];P.add_promoters = false;
G=get_refseq_introns_and_exons('hg19',P)
X=keep_fields(G,{'name','chr','tx_start','tx_end'});
flank5 = 3000;
flank3 = 0;
gidx = grep('+',G.strand,1); % plus strand
X.tx_start(gidx) = X.tx_start(gidx) - flank5;
X.tx_end(gidx) = X.tx_end(gidx) + flank3;
gidx = grepv('+',G.strand,1); % minus strand
X.tx_start(gidx) = X.tx_start(gidx) - flank3;
X.tx_end(gidx) = X.tx_end(gidx) + flank5;
X.gc = 0.50*ones(slength(X),1);
X.len = X.tx_end - X.tx_start + 1;
X.mem = ones(slength(X),1);
save_struct_noheader(X,'/xchip/cga1/lawrence/capture/hg19.full_transcript_with_3Kb_and_0Kb_flanks.20130410.txt');
% check for overlaps
X = sort_struct(X,{'chr','tx_start'});
X.gap = [nan;X.tx_end(2:end) - X.tx_start(1:end-1)];
sum(X.gap<0) % 23  (the interchromosomal jumps):   perfect!



% hg19 exons with 100bp flanks
t = load_target_file('/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.txt');
t.chr = convert_chr_back(t.chr);
t.start = t.start - 100;
t.end = t.end + 100;
t2 = condense_regions(t);
t2.gc = 0.50*ones(slength(t2),1);
t2.len = t2.end - t2.start + 1;
t2.mem = ones(slength(t2),1);
save_struct_noheader(t2,'/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_flank100bp.hg19.txt');







