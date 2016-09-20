function align_split_reads(junc,seqs,cutoff)

if ~exist('cutoff','var'), cutoff=80; end

seqs = setdiff(seqs,{'','-1'});

ns = length(seqs);
s=cell(ns,1);a=nan(1,ns);b=cell(1,ns);c=nan(2,ns);
for i=1:ns
  [aa1 bb1 cc1] = swalign(seqs{i},junc,'ALPHABET','NT');
  [aa2 bb2 cc2] = swalign(rc(seqs{i}),junc,'ALPHABET','NT');
  if aa1>aa2
    a(i)=aa1; b{i}=bb1; c(:,i)=cc1; s{i} = seqs{i};
  else
    a(i)=aa2; b{i}=bb2; c(:,i)=cc2; s{i} = rc(seqs{i});
  end
end
idx = find(a>=cutoff); [tmp ord] = sort(c(2,idx));

if length(idx)>=2
  figure(3);clf
  showalignment(multialign([junc;s(idx(ord))],'gapopen',20,'terminalGapAdjust',true));
  set(gcf,'position',[129         638        1107         309]);
end
