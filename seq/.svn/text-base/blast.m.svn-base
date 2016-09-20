function H = blast(seq, database, params, program)
%
% H = blast(seq, database, params, program)
%
% seq is sequence to blast
% database (optional) = database to search (default is 'hg18')
% params (optional) --> automatically includes '-m8'
% program (optional) = blastn (default), blastp, etc.
%
% Mike Lawrence 2008-08-27
%

if ~exist('database','var'), database = 'hg18'; end 
if ~exist('program','var'), program = 'blastn'; end
if ~exist('params','var')
  if strcmp(database,'hg18')
    if length(seq)<40
      params = '-e1e-1 -FF';
    else
      params = '-e1e-50 -nT -FF';  % nT = megablast; FF = no DUST filtering
    end
  else
    params = '';
  end
end

params = ['-m8 ' params];     % m8 = tabular form

if strcmp(database,'hg18'), database = '/xchip/tcga/gbm/analysis/lawrence/genome/hg18/orig/hg18.fa'; end

H = struct('score',[]);
try 
  [status,result] = unix(['echo ' seq ' | blastall -p ' program ' -d ' database ' ' params]);
  if status==0
    H = parse_blast_result(result);
  else
    fprintf('Problem with blast: %s\n', result);
  end
catch err
  fprintf('Problem with blast: %s\n', err.message);
end

H = sort_struct(H,'score',-1);
