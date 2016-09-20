function jobid = fishreads(sample,tn,baits,outsubdir,banner)

if ~exist('banner','var')
  banner = [sample_to_short_sample(sample) 'FISH' tn(1)];
end

try

fprintf('fishreads\n\tsample = %s\n\n',sample);
direc = ['/xchip/tcga_scratch/lawrence/' sample];

if upper(tn(1))=='T', tn = 'tumor';
elseif upper(tn(1))=='N', tn = 'normal';
else error('tn must be "t" or "n"');
end

bamfile = [direc '/' tn '.bam'];
baifile = find_bai(bamfile);
outdir = [direc '/' outsubdir];
if ~exist(outdir,'dir'), mkdir(outdir); end

baitfile = [outdir '/baits.txt'];
save_lines(baits, baitfile);

cmd = ['"java -classpath /xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
  '/xchip/tcga/gbm/analysis/lawrence/sam FishReads ' bamfile ' ' baifile ' ' baitfile ' ' outdir '"'];

jobid = bsub(cmd,banner);

catch me, excuse(me); end
