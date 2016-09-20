function F = load_fastb(file)

r = num2str(rand);
tmp_file = ['/xchip/tcga/gbm/analysis/lawrence/tmp/tmp_' r '.fasta'];
res_file = ['/xchip/tcga/gbm/analysis/lawrence/tmp/tmp_' r '.out'];

exe_file = '/wga/dev/asivache/Arachne/bin_x86_64_suse/Fastb2Fasta';
cmd = ['setenv LD_LIBRARY_PATH ' ...
       '/usr/intel/compiler81/lib:.:/util/gcc-4.3.0/lib64:' ...
       '/util/gcc-4.3.0/lib/gcc/x86_64-unknown-linux-gnu/4.3.0:/util/gcc-4.3.0/lib/:/util/lib' ...
       ';unlimit stacksize' ...
       ';' exe_file ' IN="' file '" OUT="' tmp_file '" > ' res_file];

result = system(cmd);
if result==0
  F = load_fasta(tmp_file);
else
  fprintf('\n==========================\n');
  fprintf('Fastb2Fasta had a problem!\n');
  fprintf('==========================\n');
  X = load_textfile(res_file);
  fprintf('%s',X);
  F = [];
end

% cleanup
cmd = ['rm ' tmp_file ' ' res_file];
system(cmd);
