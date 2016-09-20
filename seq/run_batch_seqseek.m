function run_batch_seqseek(samps,fasta,indir,outdir,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'queue','hour');

samps.tbam = regexprep(samps.name,'(.*)',[indir '/$1/wgs/bam/raw/tumor.bam']);
samps.nbam = regexprep(samps.name,'(.*)',[indir '/$1/wgs/bam/raw/normal.bam']);
demand_file([samps.tbam;samps.nbam]);
blacklist = '/xchip/cga2/lawrence/cga/trunk/reference/lane_blacklist.txt';
valres = 1; allowOneMismatch = 1;
for i=1:slength(samps), disp(samps.name{i});
  outstem = [outdir '/' samps.name{i}];
  keylen = getBAMFileReadLength(samps.tbam{i});
  run_seqseek(samps.tbam{i},blacklist,fasta,keylen,valres,allowOneMismatch,[outstem '-Tumor'],P);
  keylen = getBAMFileReadLength(samps.nbam{i});
  run_seqseek(samps.nbam{i},blacklist,fasta,keylen,valres,allowOneMismatch,[outstem '-Normal'],P);
end
