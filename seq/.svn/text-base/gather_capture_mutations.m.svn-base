function gather_capture_mutations(P)
% Mike Lawrence 2009-07-16

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'capture_mutation_basedir','/xchip/tcga_scratch/ng');
P=impose_default_value(P,'capture_mutation_list_output_file','*required*');

try

idx = 1;
X = {};

d = dir(P.capture_mutation_basedir);
fprintf('File: ');
for i=1:length(d)
  fname = [P.capture_mutation_basedir '/' d(i).name '/capture/mut/mutation_reports5f.maf.annotated'];
  if exist(fname,'file')
    d2 = dir(fname);
    if d2(1).bytes>0
      fprintf('%d ',idx);
      X{idx} = load_struct(fname);
      idx=idx+1;
    end
  end
end
fprintf('\n');

X = combine_structs(X);
save_struct(X,P.capture_mutation_list_output_file);

catch me; excuse(me); end
