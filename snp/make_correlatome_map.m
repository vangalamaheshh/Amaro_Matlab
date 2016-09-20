function [X,Dist]=make_correlatome_map(C,regs,fname,only_pos,pv_cutoff,dont_use_keep)

if ~exist('only_pos','var')
  only_pos=0;
end

if ~exist('dont_use_keep','var')
  dont_use_keep=0;
end

if ~exist('pv_cutoff','var')
  pv_cutoff=0.05;
end

start_supid=find_supid(C,'START');

X=reorder_D_sup(C,'cols',(start_supid+1):size(C.supdat,1));

cutoff=5;

tot=sum(X.supdat==1,2);
if ~dont_use_keep
  keep=zeros(size(X.supdat,1),1);
  keep(find(~cellfun('isempty',regexp(cellstr(X.supacc),'.*B$'))))=1;
  disp(find(~cellfun('isempty',regexp(cellstr(X.supacc),'.*B$'))));
  for k=1:2
    if (k==1) ad='AP'; else ad='DP'; end
    for i=1:length(regs{k})
      s=find_supid(X,{['Any-' ad num2str(i)],['x1-' ad num2str(i)], ['x2-' ad num2str(i)]},'exact');
      keep(s(1))=1;
      if (tot(s(2))>=cutoff) && (tot(s(3))>=cutoff) 
        keep([s(2) s(3)])=1;
      else
        keep([s(2) s(3)])=0;
      end  
    end
  end
  find(keep)
  X=reorder_D_sup(X,'cols',find(keep));
end
tot=sum(X.supdat==1,2);
X=reorder_D_sup(X,'cols',find(tot>=cutoff));

X.supdesc

nodes=[];
for i=1:size(X.supdat,1)
  nodes(i).id=['node' num2str(i)];
  tmp=deblank(regexp(X.supdesc(i,:),'([^|]*)|','tokens'));
%  [ X.supdesc(i,:) tmp{1}{1} ]
  if isempty(deblank(tmp{1}{1}))
    nodes(i).label=deblank(X.supacc(i,:));
  else
    nodes(i).label=tmp{1}{1};
  end
%  nodes(i).label=deblank(X.supacc(i,:));
end

dat=X.supdat;

Dist1=dist(dat,dat,'fisher');
Dist2=dist(dat,dat,struct('method','fisher','fdr',0.25));
Dist=Dist2;
Dist=Dist1;
Dist(abs(Dist)<-log2(pv_cutoff))=0;
Dist=Dist-diag(diag(Dist));

% remove links between x1 and x2 and Any of the same lesion
if ~isempty(regs)
  for k=1:2
    if (k==1) ad='AP'; else ad='DP'; end
    for i=1:length(regs{k})
      s=find_supid(X,{['Any-' ad num2str(i)],['x1-' ad num2str(i)], ['x2-' ad num2str(i)]},'exact');
      Dist(s,s)=0;
    end
  end
end

%Dist=Dist.^2;
Dist=Dist/max(Dist(:));
if only_pos
  Dist(Dist<0)=0;
end

write_graph_xml(fname,nodes,Dist,0);
% unix(['cp ' fname ' ~/public_html/genes/genes.fdr.xml']);
