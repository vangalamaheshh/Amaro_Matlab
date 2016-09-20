function G = convert_ttM_to_G(M)

G=[]; G.rank=(1:slength(M.gene))'; G.name = M.gene.name;

G.Nnon = sum(sum(M.Nnon(:,end,:),3),2);
G.Nsil = sum(sum(M.Nsil(:,end,:),3),2);
G.Nflank = sum(sum(M.Nflank(:,end,:),3),2);

G.nnon = sum(sum(M.nnon,3),2);
G.nsil = sum(sum(M.nsil,3),2);
G.nflank = sum(sum(M.nflank,3),2);
