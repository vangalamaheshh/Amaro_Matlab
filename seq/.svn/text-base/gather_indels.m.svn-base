function gather_indels(P)
% Mike Lawrence 2009-07-16

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'indel_basedir','/xchip/tcga_scratch/ng');
P=impose_default_value(P,'nextdir_mask','*');
P=impose_default_value(P,'indel_subdir','capture/indel');
P=impose_default_value(P,'indel_list_output_file','*required*');

try

idx = 1;
X = {};

d = dir([P.indel_basedir '/' P.nextdir_mask]);
keyboard
fprintf('File: ');
for i=1:length(d)
  d2 = [P.indel_basedir '/' d(i).name '/' P.indel_subdir];
  fd = dir([d2 '/*.somatic.txt']);
  if ~isempty(fd)
    if length(fd)>1, error('Confused: multiple possible files in %s',d2); end
    fname = [d2 '/' fd.name];
    d2 = dir(fname);
    if d2(1).bytes>0
      fprintf('%d ',idx);
      X{idx} = load_indel_file(fname);
      X{idx}.file = repmat({fname},slength(X{idx}),1);
      idx=idx+1;
    end
  end
end
fprintf('\n');

X = combine_structs(X);
X = keep_fields(X,{'file','chr','start','end','details1','status','type','gene'});

save_struct(X,P.indel_list_output_file,'no_headers');

catch me; excuse(me); end
