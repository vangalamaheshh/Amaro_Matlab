function create_blast_db(fasta_fname,dirname,dbname,params)
% create_blast_db(fasta_fname,dirname)
%
% creates a blast database named <dbname> in <dirname>,
% using the sequences in fasta file <fasta_fname>
% optional: can specify additional params in <params>
%
% Mike Lawrence 2009-03-03

if ~exist('fasta_fname','var'), error('Must specify input fasta filename'); end
if ~exist('dirname','var'), error('Must specify output directory'); end
if ~exist('dbname','var'), error('Must specify output database name'); end
if ~exist('params','var'), params=[]; end

exe_name = '/xchip/tcga/gbm/analysis/lawrence/blast/bin/formatdb';
params = ['-p F ' params];

cmd = [exe_name ' -i ' fasta_fname ' ' params];
[result output] = system(cmd);

if result~=0, error('Error creating blast database: %s',result); end

