function [in_12,in_1,in_2]=textdiff(t1,t2)

[Mt,m1,m2]=match_string_sets_hash(t1,t2);

in_1=t1(setdiff(1:length(t1),m1));
in_2=t2(setdiff(1:length(t2),m2));
in_12=t1(unique(m1));

