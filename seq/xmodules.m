function [keep_genes keep_pats m v x e labs stats] = xmodules(M,outfname,silent)
if ~exist('outfname','var') || isempty(outfname)
  outfile=1;
else
  outfile=fopen(outfname,'w');
end

if ~exist('silent','var')
  silent=0;
end

require_fields(M,{'gene','patient'});

[gene gi gj] = unique(M.gene);
[patient pi pj] = unique(M.patient);

S=sparse(pj,gj,ones(length(pj),1));

muts_per_gene = histc(gj,1:length(gene));
keep_genes = find(sum(S,1)>1);
keep_pats = find(sum(S,2)>1);

S2=S(keep_pats,keep_genes);
[pj2 gj2]=find(S2);
m=[pj2 gj2];

v=[];
for i=1:length(keep_pats)
  idx=find(m(:,1)==i);
  for c1=1:length(idx)-1
    for c2=(c1+1):length(idx)
      v=[v; idx(c1) idx(c2)];
    end
  end
end

gv=[m(v(:,1),2) m(v(:,2),2)];
[ugv,ui,uj]=unique(gv,'rows');

x=[];
for i=1:length(ui)
  idx=find(uj==i);
  for c1=1:length(idx)-1
    for c2=(c1+1):length(idx)
      x=[x; v(idx(c1),:) v(idx(c2),:)];
    end
  end
end

e=[];
[um,umi,umj]=unique(x(:));
for i=1:length(umi)
  idx=find(umj==i);
  for j=1:length(idx)-1
    e=[e; mod(idx(j)-1,size(x,1))+1 mod(idx(j+1)-1,size(x,1))+1 ];
  end
end

labs=unionfind(e,size(x,1));

stats=[];
for i=1:max(labs)
  idx=find(labs==i);
%  if length(idx)<2, continue; end;
  x_in_mod=x(idx,:);
  g_in_mod = keep_genes(unique(m(unique(x_in_mod(:)),2)));
  p_in_mod = keep_pats(unique(m(unique(x_in_mod(:)),1)));
  stats.num_x(i)=length(idx);
  stats.num_mut(i)=length(unique(x_in_mod(:)));
  stats.num_genes(i)=length(g_in_mod);
  stats.num_pats(i)=length(p_in_mod);
%  disp(i);
%  disp(gene(g_in_mod));
%  disp(patient(p_in_mod));
%  disp(' ');
  if ~silent
    mS=S(p_in_mod,g_in_mod);
%    tot_r=sum(mS,2); [tmp, sr]=sort(tot_r,'descend');
%    tot_c=sum(mS,1); [tmp, sc]=sort(tot_c,'descend');
%    mS=mS(sr,sc);
    for ri=0:length(p_in_mod)
      for ci=0:length(g_in_mod)
        if (ri==0) && (ci==0)
          fprintf(outfile,'%s\t',['Module #' num2str(i)]); 
        elseif ri==0
          fprintf(outfile,'%s\t',gene{g_in_mod(ci)});
        elseif ci==0
          fprintf(outfile,'%s\t',patient{p_in_mod(ri)});
        elseif mS(ri,ci)
          fprintf(outfile,'x\t');
        else
          fprintf(outfile,'\t');
        end          
      end
      fprintf(outfile,'\n');
    end
    fprintf(outfile,'\n\n');
  end
end

if outfile~=1
  fclose(outfile);
end

