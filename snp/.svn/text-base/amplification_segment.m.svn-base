function tbl=amplification_segment(C,rg,cyto,refseq_genes,ts,dt)
%
%---
% $Id$
% $Date: 2007-12-19 09:49:46 -0500 (Wed, 19 Dec 2007) $
% $LastChangedBy: rameen $
% $Rev$

nms={rg(:).refseq};
[M,m1,m2]=match_string_sets_hash(nms,refseq_genes);

tbl=cell(length(refseq_genes),size(C.dat,2));

for i=1:length(m1)
  if mod(i,100)==0
    disp(i);
  end
  p{i}=find_snps(C,chromosome2num(rg(m1(i)).chr(4:end)),rg(m1(i)).start,rg(m1(i)).end,1);
  chr=chromosome2num(rg(m1(i)).chr(4:end));
  if ~isnan(chr)
    in_chr=find(C.chrn==chr);
  else
    continue;
  end
  if isempty(in_chr)
    continue;
  end
  pc{i}=p{i}-min(in_chr)+1;
  x=C.dat(in_chr,:);
  if ~isempty(p{i})
    vals=mean(C.dat(p{i},:),1);
    x=x.*repmat(sign(vals),size(x,1),1);
    tv=ts(floor(1.5-0.5*sign(vals))); % ts=[amp_t del_t]
    dtv=dt(floor(1.5-0.5*sign(vals))); % dt=[amp_dt del_dt]
    vals=abs(vals);
    tv=abs(tv);
    rls=runlength(x>repmat(max(vals-dtv,tv),size(x,1),1));
    for j=1:size(C.dat,2)
      if vals(j)>tv(j)
        seg_i=find(rls{j}(:,1)<=max(pc{i}) & rls{j}(:,2)>=min(pc{i}));
        if rls{j}(seg_i,3)~=1
          disp('BAD');
          keyboard
        end
        tbl{i,j}=genomic_location(C,{in_chr(rls{j}(seg_i,1):rls{j}(seg_i,2))},cyto,1);
      end
    end
  end  
end

