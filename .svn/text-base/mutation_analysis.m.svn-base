function gene_rest=mutation_analysis(D,gene_supacc,do_perm)
cd ~/projects/snp
genem_supid=strmatch(gene_supacc,D.supacc,'exact');
gene_rest.cls0=find(D.supdat(genem_supid,:)==0);
gene_rest.cls1=find(D.supdat(genem_supid,:)==1);
if ~exist('do_perm','var')
  do_perm=1;
end
[gene_rest.ap2,gene_rest.p2,gene_rest.p1_gte,gene_rest.p1_lte,gene_rest.sets,...
        gene_rest.s,gene_rest.n,gene_rest.acc,gene_rest.select_compounds1,gene_rest.select_compounds2]=...
        generate_data_for_figure(D,gene_rest.cls0,gene_rest.cls1,gene_supacc,do_perm);

figure(1); clf;
subplot(2,1,1);
%imagesc(D.dat([gene_rest.select_compounds1 gene_rest.select_compounds2],[gene_rest.cls0 gene_rest.cls1]));
bluepink;
imagesc_trim(-dna_norm(D.dat([gene_rest.select_compounds1 gene_rest.select_compounds2],[gene_rest.cls0 gene_rest.cls1])));
bluepink;
subplot(2,1,2);
image((D.supdat(5:9,[gene_rest.cls0 gene_rest.cls1])+0)*32+32);
set(gca,'YTick',1:size(D.supdat,1));
set(gca,'YTickLabel',D.supacc(5:9,:));
D.gacc(gene_rest.select_compounds1,:)
keyboard
print('-dpng','-f1','-r180',[gene_supacc '_comp.png']);

x=gene_rest;
f=fopen([gene_supacc '_comp.txt'],'w');
for i=1:length(x.select_compounds1)
  fprintf(f,'%s\t%e\t%f\r\n',D.gacc(x.select_compounds1(i),:),x.ap2(x.select_compounds1(i)), ...
          x.s(x.select_compounds1(i)));
end
for i=length(x.select_compounds2):-1:1
  fprintf(f,'%s\t%e\t%f\r\n',D.gacc(x.select_compounds2(i),:),x.ap2(x.select_compounds2(i)), ...
          x.s(x.select_compounds2(i)));
end
fclose(f);


