function I = find_near_introns(E,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'range_begin',10);
P = impose_default_value(P,'range_end',100);

[g gi gj] = unique(E.gene);
ng = length(g);

fprintf('Finding near introns: ');
I = cell(ng,1);
for i=1:ng, if ~mod(i,1000), fprintf('%d/%d ',i,ng); end
  idx = find(gj==i);
  p = [E.chr(idx) E.start(idx) E.end(idx)];
  [c ci cj] = unique(p(:,1));
  allin = [];
  for chri=1:length(c)
    jdx = find(cj==chri);
    ex = sortrows(p(jdx,2:3));
    nex = size(ex,1);
    inl = [ex(:,1)-P.range_end ex(:,1)-P.range_begin];
    inr = [ex(:,2)+P.range_begin ex(:,2)+P.range_end];
    for j=1:nex
      if j>1
        inl(j,1) = max(inl(j,1),inr(j-1,2)+1);
        inl(j,2) = max(inl(j,2),inr(j-1,2)+1);
      end
      if j<nex
        inr(j,1) = min(inr(j,1),ex(j+1,1)-P.range_begin);
        inr(j,2) = min(inr(j,2),ex(j+1,1)-P.range_begin);
      end
    end
    in = sortrows([inl; inr]);
    in = in(in(:,2)>in(:,1),:);
    in = [c(chri)*ones(size(in,1),1) in];
  end
  allin = [allin; in];
  I{i}.gene = repmat({g{i}},size(allin,1),1);
  I{i}.chr = allin(:,1);
  I{i}.start = allin(:,2);
  I{i}.end = allin(:,3);
end, fprintf('\n');

I = concat_structs(I);

I.gc = get_gc_content(I.chr,I.start,I.end,P.build);
I.len = I.end-I.start+1;
