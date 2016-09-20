function [lnk,idx]=blank_lnk_idx(n)

idx=1:n;
lnk=[ ones(n-1,1) (1:(n-1))' (2:n)' (2:n)' zeros(n-1,1)];
