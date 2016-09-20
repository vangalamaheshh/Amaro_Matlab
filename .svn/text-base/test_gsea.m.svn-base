n=10000;

r1=randperm(n);

np=1000;
sl=50;
for i=1:np
  r=randperm(n);
  sets{i}=r(1:sl);
end

[ks,kspv,rs,rspv]=gsea(r1,sets);

[tmp,revorder]=sort(r1);
figure(1);clf;

plot_gsea(revorder(sets{607}),n);

dat=rand(n,50);
sup=[ ones(10,1); zeros(40,1)];

[trueres,randres,permpv]=gsea_phen_perm(dat,sup,'ttest',sets(1:100),20);

sets_path='~/projects/snp/gene_sets/';
go_sets=read_mit_gmt_file([ sets_path 'HG_U95Av2_GO.gmt' ]);
biocarta_sets=read_mit_gmx_file([sets_path 'HG_U95Av2_biocarta.gmx']);
cyto_sets=read_mit_gmt_file([ sets_path 'HG_U95Av2_cyto.gmt' ]);
gearray_sets=read_mit_gmx_file([sets_path 'HG_U95Av2_gearray.gmx']);
genmapp_sets=read_mit_gmt_file([ sets_path 'HG_U95Av2_genmapp.gmt' ]);
proteome_sets=read_mit_gmt_file([ sets_path 'HG_U95Av2_proteome.gmt' ]);

all_sets=[ go_sets biocarta_sets cyto_sets gearray_sets genmapp_sets ...
           proteome_sets];

save([sets_path 'all_sets.mat'],'all_sets','go_sets','biocarta_sets',...
      'cyto_sets','gearray_sets','genmapp_sets','proteome_sets');








