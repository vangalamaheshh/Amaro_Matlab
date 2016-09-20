function d=dna_norm_con(v,con_cls_array,make_like_cls)

d=zeros(size(v));
for i=1:length(con_cls_array)
  [d(:,con_cls_array{i}),m{i},s{i}]=dna_norm_nan(v(:,con_cls_array{i}));
end
if nargin>2
  d=d.*repmat(s{make_like_cls},1,size(d,2))+ ...
    repmat(m{make_like_cls},1,size(d,2));
end

