function [table pvalue]=construct_fishers(X,field1,field2)

a=sum(X.(field1) & X.(field2));
b=sum(X.(field1) & ~X.(field2));
c=sum(~X.(field1) & X.(field2));
d=sum(~X.(field1) & ~X.(field2));

%pvalue=myfisher([a,b;c,d]);
table=[a,b;c,d];
[pvalue,p2]=fisher_exact_test(a,b,c,d,1e-20,'right');

end






function test

cd /Users/stewart/Projects/Cancer/mouse
load('synteny_workspace.mat')
% PFDEL=load_struct_noheader('/Users/amaro/Downloads/PfiferGenesDel.txt');
% PFAMP=load_struct_noheader('/Users/amaro/Downloads/PfiferGenesAmp.txt');
% PFAMP_genes=unique(PFAMP.col13);
% PFDEL_genes=unique(PFDEL.col13);
% gene_order.PfAmp=ismember(gene_order.human,PFAMP_genes);
% gene_order.PfDel=ismember(gene_order.human,PFDEL_genes);
% MouseFinalDeletions=load_struct('/Users/amaro/Documents/uniqMouseDeletions.txt');
% MouseFinalAmplifications=load_struct('/Users/amaro/Documents/uniqMouseAmplifications.txt');
% gene_order.MouseAmp=ismember(gene_order.mouse,MouseFinalAmplifications.MouseAmplifications);
% gene_order.MouseDel=ismember(gene_order.mouse,MouseFinalDeletions.MouseDeletions);
% 
% refGene=load_struct_noheader('/Users/amaro/Downloads/refGeneFile.txt');
% refGeneMouse=load_struct_noheader('/Users/amaro/Downloads/refGeneFile_Mouse.txt');
% for i=1:slength(gene_order)
% loc=find(ismember(refGene.col13,gene_order.human{i}),1,'first');
% gene_order.HumanCoord{i,1}=sprintf('%s-%s',refGene.col3{loc},refGene.col5{loc});
% loc=find(ismember(refGeneMouse.col13,gene_order.mouse{i}),1,'first');
% gene_order.MouseCoord{i,1}=sprintf('%s-%s',refGeneMouse.col3{loc},refGene.col5{loc});
% if mod(i,1000)==0
% i
% end
% end



OutFileAmp.CommonlyAmplified=gene_order.human(gene_order.PfAmp);
OutFileAmp.HumanCoordAmp=gene_order.HumanCoord(gene_order.PfAmp);
OutFileAmp.MouseAmp=gene_order.mouse(gene_order.PfAmp);
OutFileAmp.MouseCoordAmp=gene_order.MouseCoord(gene_order.PfAmp);
OutFileAmp.InMouseAmp=(gene_order.MouseAmp(gene_order.PfAmp));
save_struct(OutFileAmp,'PfiferAmps.txt');
OutFileDel.CommonlyDeleted=gene_order.human(gene_order.PfDel);
OutFileDel.HumanCoordDel=gene_order.HumanCoord(gene_order.PfDel);
OutFileDel.MouseDel=gene_order.mouse(gene_order.PfDel);
OutFileDel.MouseCoordDel=gene_order.MouseCoord(gene_order.PfDel);
OutFileDel.InMouseDel=(gene_order.MouseDel(gene_order.PfDel));
save_struct(OutFileDel,'PfiferDel.txt');




OutFileAmp.CommonlyAmplified=gene_order.human(gene_order.RudinA);
OutFileAmp.HumanCoordAmp=gene_order.HumanCoord(gene_order.RudinA);
OutFileAmp.MouseAmp=gene_order.mouse(gene_order.RudinA);
OutFileAmp.MouseCoordAmp=gene_order.MouseCoord(gene_order.RudinA);
OutFileAmp.InMouseAmp=(gene_order.MouseAmp(gene_order.RudinA));
save_struct(OutFileAmp,'RudinAmps.txt');
OutFileDel.CommonlyDeleted=gene_order.human(gene_order.RudinD);
OutFileDel.HumanCoordDel=gene_order.HumanCoord(gene_order.RudinD);
OutFileDel.MouseDel=gene_order.mouse(gene_order.RudinD);
OutFileDel.MouseCoordDel=gene_order.MouseCoord(gene_order.RudinD);
OutFileDel.InMouseDel=(gene_order.MouseDel(gene_order.RudinD));
save_struct(OutFileDel,'RudinDel.txt');






gene_order

[t1,p(1)]=construct_fishers(gene_order,'MouseDel','RudinD')
[t2,p(2)]=construct_fishers(gene_order,'MouseAmp','RudinA')
[t3,p(3)]=construct_fishers(gene_order,'MouseDel','PfDel')
[t4,p(4)]=construct_fishers(gene_order,'MouseAmp','PfAmp')

fprintf('%.3e\n',p)

char(unique(gene_order.human(gene_order.MouseDel)))
char(unique(gene_order.human(gene_order.RudinD)))
char(unique(gene_order.human(gene_order.PfDel)))

char(unique(gene_order.human(gene_order.MouseAmp)))
char(unique(gene_order.human(gene_order.RudinA)))
char(unique(gene_order.human(gene_order.PfAmp)))

% -> venny 

end