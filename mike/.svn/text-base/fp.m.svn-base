function fp(adj)
% fix plot

if ~exist('adj','var'), adj=0.05; end

pos = get(gca,'position');
pos(2)=pos(2) + adj;
pos(4)=pos(4) - 0.5*adj;
set(gca,'position',pos);
