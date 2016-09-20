
N=10;
d=rand(N,N);
D.gdesc=num2str((1:N)');
D.gacc=D.gdesc;
D.sdesc=D.gdesc;
D.dat=[ ones(N/2,N/2) zeros(N/2,N/2);zeros(N/2,N/2) ones(N/2,N/2)]+d;
r1=randperm(N);
r2=randperm(N);
D1=reorder_D_rows(reorder_D_cols(D,r1),r2);
write_mit_gct_file('demo.gct',D1);

gp = GenePatternServer('http://genepattertest:8080','gadgetz@broad.mit.edu');

res=runAnalysis(gp,

gp_clustering_figure -o demo.reord.reord.gct.dendr -d demo.reord.reord.gct -dp 0 -p -c -cc -gd demo.dend.odf -gdr 0 -sd demo.reord.dend.odf -sdr 0 -ot 0 -sf 5 -gf 10 -sx 0.4 0 1.2 -sy 0.2 0.05 0 0.4 0.05 -l 1 -lnc 2 -ltf 12 -lf 6
