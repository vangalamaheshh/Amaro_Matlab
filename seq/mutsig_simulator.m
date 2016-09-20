function G = mutsig_simulator(G,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'genes_to_analyze',1:slength(G));

ng = slength(G);
z = nan(ng,1); G.nnon=z; G.nsil=z; G.nflank=z;
for i=1:ng, modi(i,100);
  if ~ismember(i,P.genes_to_analyze), continue; end


  G.nnon(i) = round(G.Nnon(i)*G.ground_truth(i));
  G.nsil(i) = round(G.Nsil(i)*G.ground_truth(i));
  G.nflank(i) = round(G.Nflank(i)*G.ground_truth(i));


%  G.nnon(i) = binornd(G.Nnon(i),G.ground_truth(i));
%  G.nsil(i) = binornd(G.Nsil(i),G.ground_truth(i));
%  G.nflank(i) = binornd(G.Nflank(i),G.ground_truth(i));


%  G.nnon(i) = binornd_fast(G.Nnon(i),G.ground_truth(i));
%  G.nsil(i) = binornd_fast(G.Nsil(i),G.ground_truth(i));
%  G.nflank(i) = binornd_fast(G.Nflank(i),G.ground_truth(i));
end

