function F = extract_from_fastb(file,ids)
% extract_from_fastb(file,ids)
%
% calls SelectFastb to extract the specified sequences from the specified fastb file.
% ids are the zero-based indices of the desired sequences in the file
%
% returns a fasta structure
%
% Mike Lawrence 2008-11-18

in_file = file;
r = num2str(rand);
id_file = ['/xchip/tcga/gbm/analysis/lawrence/tmp/tmp_id_' r '.txt'];
out_file = ['/xchip/tcga/gbm/analysis/lawrence/tmp/tmp_out_' r '.txt'];
exe_file = '/wga/dev/asivache/Arachne/bin_x86_64_suse/SelectFastb';

[ids ids_i ids_j] = unique(ids);
% HAVE TO DO THIS BECAUSE OF A BUG IN SelectFastB
%    which causes it to crash in certain cases when ids are duplicated

out=fopen(id_file,'wt');
for i=1:length(ids), fprintf(out,'%d\n',ids(i)); end
fclose(out);

cmd = ['setenv LD_LIBRARY_PATH ' ...
       '/usr/intel/compiler81/lib:.:/util/gcc-4.3.0/lib64:' ...
       '/util/gcc-4.3.0/lib/gcc/x86_64-unknown-linux-gnu/4.3.0:/util/gcc-4.3.0/lib/:/util/lib' ...
       ';unlimit stacksize' ...
       ';' exe_file ' PRE="" IN="' in_file '" IDS="' id_file '" > ' out_file];

try
  result = system(cmd);
catch
  result = -1;
end

if result==0
  F = load_fasta(out_file);
  F = reorder_struct(F,ids_j);
else
  fprintf('\n==========================\n');
  fprintf('SelectFastb had a problem!\n');
  fprintf('==========================\n');
  x=load_textfile(out_file);
  fprintf('%s',x);
  F = [];
end

% cleanup

cmd = ['rm ' id_file ' ' out_file];

try
  system(cmd);
catch
  fprintf('Unable to delete file\n');
end

