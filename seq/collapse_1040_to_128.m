function out = collapse_1040_to_128(in)
% input = 1040 categories from "allcateg"
%
% output = 64 untranscribed + 64 transcribed

if ~isnumeric(in), error('input should be numeric'); end

c1040 = load_struct('/xchip/tcga_scratch/lawrence/db/allcateg/categs.txt');
c128 = load_struct('/xchip/tcga_scratch/lawrence/db/allcateg/categs128.txt');

tmp1040 = parse(c1040.name,'^(.*):(.*):(.*)$',{'tx','zone','ctxt'});
tmp128 = parse(c128.name,'^(.*):(.*)$',{'tx','ctxt'});

idx_nontx = find(strcmp('neither_transcribed',tmp1040.tx));
idx_tx = setdiff(1:1040,idx_nontx);

out = nan(128,size(in,2));
for i=1:128
  idx1 = find(strcmp(tmp1040.ctxt,tmp128.ctxt{i}));
  switch tmp128.tx{i}
    case 'nontranscribed'
      idx = intersect(idx1,idx_nontx);
    case 'transcribed'
      idx = intersect(idx1,idx_tx);
  end
  out(i,:) = sum(in(idx,:));
end

