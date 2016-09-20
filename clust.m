function [y,Dord]=clust(x)

D=make_D(x);
Dord=two_way_clustering(D);
y=Dord.dat;
