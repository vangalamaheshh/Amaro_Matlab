function jobs = run_seqseek(bam,blacklist,fasta,keylen,valres,allowOneMismatch,outstem,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'queue','hour');

jcp = get_jcp;
class = 'org.broadinstitute.cga.tools.seq.SeqSeek';

cmds={};
banners={};

for c=1:24
  outfile1 = [outstem '.chr' num2str(c) '.run.txt'];
  outfile2 = [outstem '.chr' num2str(c) '.weird.txt'];
  if ~exist(outfile1,'file') || ~exist(outfile2,'file')
    cmds{end+1} = ['java -Xmx2g -classpath ' jcp ' ' class ' ' bam ' ' blacklist ' ' fasta ' '...
                 num2str(keylen) ' ' num2str(valres) ' ' num2str(allowOneMismatch) ' ' num2str(c) ' ' outstem];
    banners{end+1} = [bam '-' num2str(c)];
  end
end

%jobs = bsub(cmds,banners);
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
batchsize=6;
jobs = [];
for i=1:batchsize:length(cmds)
  idx = i:min(length(cmds),i+batchsize-1);
  jobs = [jobs;bsub(cmds(idx),banners(idx),P)];
end
