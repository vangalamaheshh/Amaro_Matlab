function single_screen_numbers(fname,r0,mut_rate,fnr,n,Ng,rnk,q,L)

if isempty(fname)
  f=1;
else
  f=fopen(fname,'w');
end

fprintf(f,'Background model\n');
fprintf(f,'========================\n');
fprintf(f,['mu=' num2str(r0/L) newline]);
fprintf(f,['average gene length=' num2str(L) newline]);
x=0:20;
fprintf(f,'-------------------------\n');
fprintf(f,'Screen\n');
fprintf(f,'-------------------------\n');
fprintf(f,['N=' num2str(n) newline]);
p=binopdf(x,n*L,r0/L);

fprintf(f,['Sequencing ' num2str(Ng) ' genes' newline]);
fprintf(f,['Prob. to have >0 mutations in screen: ' num2str(1-p(1)) ' (~' num2str(1-exp(-r0*n)) ')' newline]);
fprintf(f,['Expected number of genes with >0 mutations in discovery screen: ' num2str((1-p(1))*Ng) ' (~' num2str((1- ...
                                                  exp(-r0*n))*Ng) ')' newline]);

fprintf(f,'-------------------------\n');
fprintf(f,'Total Cost\n');
fprintf(f,'-------------------------\n');
fprintf(f,'Number of exons per gene=15\n');
fprintf(f,'Cost per exon=$2\n');
fprintf(f,['Total genes sequenced=' num2str(Ng*n) newline]); 
fprintf(f,['Total cost=' num2str((Ng*n)*30) newline]);

fprintf(f,'-------------------------\n');
fprintf(f,'Power Calculation\n');
fprintf(f,'-------------------------\n');
fprintf(f,['False negative rate=' num2str(fnr) newline]);
F=1-fnr;
fprintf(f,['Aim to detect mutations that are in ' num2str(mut_rate) newline]);
fprintf(f,['Taking false negative rate we want to detect mutations that occur at ' num2str(mut_rate*F) newline]);
r1=r0+(mut_rate*F);
fprintf(f,['Assuming ' num2str(rnk-1) ' genes are already known to be mutated' newline]);
fprintf(f,['A gene is called significant if it''s p-value is below ' num2str(q) 'x' num2str(rnk) '/' num2str(Ng) '=' ...
      num2str(q*rnk/Ng) newline]);

for i=1:length(r1)
  [po(i),ks(i)]=single_screen_power(r0,r1(i),n,Ng,rnk,q,L);
  fprintf(f,['Power to detect ' num2str(mut_rate(i)) ' is ' num2str(po(i)) ' (critical x=' num2str(ks(i)) ')' newline]);
end

if f~=1
  fclose(f);
end
