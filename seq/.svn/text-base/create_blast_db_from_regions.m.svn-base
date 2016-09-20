function create_blast_db_from_regions(R,dirname,dbname,params,build)
% create_blast_db_from_regions(R,dirname,dbname,params,build)
%
% creates a blast database named <dbname> in <dirname>,
% using the sequences specified in struct R (with fields chr,start,end
% optional: can specify additional params in <params>
% build = hg18 by default
%
% Mike Lawrence 2009-03-03

if ~exist('R','var'), error('Must specify input regions'); end
if ~exist('dirname','var'), error('Must specify output directory'); end
if ~exist('dbname','var'), error('Must specify output database name'); end
if ~exist('params','var'), params=[]; end
if ~exist('build','var'), build='hg18'; end

require_fields(R,{'chr','start','end'});

F=[];
for i=1:slength(R)
  c = R.chr(i);
  s = R.start(i);
  e = R.end(i);
  if iscell(c), c=c{1}; end
  if isnumeric(c), F.header{i,1} = sprintf('chr%d:%d-%d',c,s,e);
  else F.header{i,1} = sprintf('%s:%d-%d',c,s,e); end
  F.seq{i,1} = genome_region(c,s,e,build);
end

if ~exist(dirname,'dir'), mkdir(dirname); end
fasta_fname = [dirname '/' dbname];
save_fasta(F,fasta_fname);
create_blast_db(fasta_fname,dirname,dbname,params);
