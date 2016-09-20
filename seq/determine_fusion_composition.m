function [desc details] = determine_fusion_composition(seq,chr1,pos1,chr2,pos2,P)

if ~exist('P','var'), P=[];end

P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'genome_radius_to_examine',400);
P = impose_default_value(P,'swalign_scoringmatrix','NUC1');
P = impose_default_value(P,'min_alignment_length_for_match',30);

seq = upper(seq);

if any(seq=='N')
%  fprintf('can''t deal with N''s yet... replacing with A''s\n');
  seq(seq=='N')='A';
end



g1 = upper(genome_region(chr1,pos1-P.genome_radius_to_examine,pos1+P.genome_radius_to_examine,P.build));
g2 = upper(genome_region(chr2,pos2-P.genome_radius_to_examine,pos2+P.genome_radius_to_examine,P.build));

if any(g1=='N')
  g1(g1=='N')='A';
end
if any(g2=='N')
  g2(g2=='N')='A';
end


swparams = {'gapopen',10,'extendgap',10,'alphabet','nt','scoringmatrix',P.swalign_scoringmatrix};
f=[]; [f.a f.b f.c] = swalign(seq,g1,swparams{:});
r=[]; [r.a r.b r.c] = swalign(seq,rc(g1),swparams{:});
if f.a>r.a, end1=f; end1.str=0; else end1=r; end1.str=1; end
end1.seqlen = sum(end1.b(1,:)~='-');
end1.reflen = sum(end1.b(3,:)~='-');
f=[]; [f.a f.b f.c] = swalign(seq,g2,swparams{:});
r=[]; [r.a r.b r.c] = swalign(seq,rc(g2),swparams{:});
if f.a>r.a, end2=f; end2.str=0; else end2=r; end2.str=1; end
end2.seqlen = sum(end2.b(1,:)~='-');
end2.reflen = sum(end2.b(3,:)~='-');

end1.chr = decell(convert_chr_back(chr1));
end1.pos = pos1;
end2.chr = decell(convert_chr_back(chr2));
end2.pos = pos2;

if end1.str==0
  end1.refstart = end1.pos-1+end1.c(2)-P.genome_radius_to_examine;
  end1.refend = end1.refstart + end1.reflen - 1;
else
  end1.refend = end1.pos+P.genome_radius_to_examine-end1.c(2)+1;
  end1.refstart = end1.refend - end1.reflen + 1;
end
if end2.str==0
  end2.refstart = end2.pos-1+end2.c(2)-P.genome_radius_to_examine;
  end2.refend = end2.refstart + end2.reflen - 1;
else
  end2.refend = end2.pos+P.genome_radius_to_examine-end2.c(2)+1;
  end2.refstart = end2.refend - end2.reflen + 1;
end

if end1.c(1)<end2.c(1), lft=end1; lft.end=1; rgt=end2; rgt.end=2;
else lft=end2; lft.end=2; rgt=end1; rgt.end=1; end

details = [];
details.end1 = end1;
details.end2 = end2;
details.lft = lft;
details.rgt = rgt;

% see if it's successful
lenok = (end1.seqlen>=P.min_alignment_length_for_match & end2.seqlen>=P.min_alignment_length_for_match);

if lenok
  details.success = true;
  % summarize as text description
  symb ={'(+)','(-)'}';
  lft.desc = ['[' lft.chr symb{lft.str+1} num2str(lft.refstart) ':' num2str(lft.refend) ']'];
  lft.prefix = lft.c(1)-1;
  if lft.prefix>0, lft.desc = ['[' num2str(lft.prefix) ']' lft.desc]; end
  rgt.desc = ['[' rgt.chr symb{rgt.str+1} num2str(rgt.refstart) ':' num2str(rgt.refend) ']'];
  rgt.suffix = length(seq)-(rgt.c(1)+rgt.seqlen)+1;
  if rgt.suffix>0, rgt.desc = [rgt.desc '[' num2str(rgt.suffix) ']']; end
  details.spacer = length(seq)-lft.seqlen-rgt.seqlen-lft.prefix-rgt.suffix;
  desc = [lft.desc '[' num2str(details.spacer) ']' rgt.desc];
else
  details.success = false;
  desc = 'alignment_failed';
  details.spacer = nan;
end

details.desc = desc;
