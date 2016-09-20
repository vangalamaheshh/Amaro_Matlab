function [v,vlst]=venn_values(lst1,lst2)

both12=intersect(lst1,lst2);
only1=setdiff(lst1,both12);
only2=setdiff(lst2,both12);

v=[length(only1) length(both12) length(only2)];
vlst={only1,both12,only2};
