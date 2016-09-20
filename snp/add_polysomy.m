function [D,supids]=add_polysomy(D,cyto,t,p)

if ~exist('p','var')
  p=0.5;
end

if ~exist('t','var')
  t=0.3;
end

if ~isfield(D,'chrn')
  D=add_chrn(D);
end

if ischar(cyto)
  cyto=read_cytoband_data(cyto);
end
nms=strvcat(cyto(:).name);

supids=[];
% full chromosomes
for is_del=0:1
  for i=1:max(D.chrn)
    pos=find(cat(1,cyto(:).chrn)==i);
    st=cyto(min(pos)).start;
    en=cyto(max(pos)).end;
    pos_snps=find_snps(D,i,st,en,1);
    if ~isempty(pos_snps)
      if ~is_del
        v=mean(D.dat(pos_snps,:)>t,1)>p;
        [D,si]=add_D_sup(D,[num2str(i) 'A'],[num2str(i) 'A'],v,'cols');
      else
        v=mean(D.dat(pos_snps,:)<-t,1)>p;
        [D,si]=add_D_sup(D,[num2str(i) 'D'],[num2str(i) 'D'],v,'cols');
      end
      supids=[supids; si];
    else
      disp(['no SNPs in ' num2str(i) ]);
    end
  end
end

% chromosome arms
for is_del=0:1
  for i=1:max(D.chrn)
    for j=['p' 'q']
      pos=find(cat(1,cyto(:).chrn)==i & nms(:,length(num2chromosome(i))+1)==j);
      st=cyto(min(pos)).start;
      en=cyto(max(pos)).end;
      pos_snps=find_snps(D,i,st,en,1);
%      [ is_del i st en double(j) length(pos_snps) ]
      if ~isempty(pos_snps)
        if ~is_del
          v=mean(D.dat(pos_snps,:)>t,1)>p;
          [D,si]=add_D_sup(D,[num2str(i) j 'A'],[num2str(i) j 'A'],v,'cols');
        else
          v=mean(D.dat(pos_snps,:)<-t,1)>p;
          [D,si]=add_D_sup(D,[num2str(i) j 'D'],[num2str(i) j 'D'],v,'cols');
        end
        supids=[supids; si];
      else
        disp(['no SNPs in ' num2str(i) j]);
      end
    end
  end
end
