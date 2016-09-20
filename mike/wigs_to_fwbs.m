function wigs_to_fwbs(indir,outdir,P)
% convert all WIGs to FWBs

if ~exist('P','var'), P = []; end
P = impose_default_value(P,'batchsize',1);   % how many conversions per LSF job
P = impose_default_value(P,'queue','hour');

ede(outdir);

tmp = pwd; cd(indir); d = direc('*wig*'); cd(tmp);
cmds = {}; jcp = get_jcp;
jcn = 'org.broadinstitute.cga.tools.seq.Wig2Fwb';
for i=1:length(d)
  infile = [indir '/' d{i}];
  outfile = [outdir '/' d{i} '.fwb'];
  if exist(outfile,'file'), continue; end
  cmds{end+1,1} = ['java -classpath ' jcp ' ' jcn ' ' infile ' 1 ' outfile];
end

nc = length(cmds); nb = ceil(nc/P.batchsize);bcmds = cell(nb,1); bbanners = cell(nb,1);
for i=1:nb
  cmdfile = [outdir '/cmds.' num2str(i) '.bsub'];
  save_lines(cmds((i-1)*P.batchsize+1:min(nc,i*P.batchsize)),cmdfile);
  bcmds{i} = ['source ' cmdfile]; bbanners{i} = ['CP' num2str(i)];
end

bsub(bcmds,bbanners,P);
