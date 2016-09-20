function [nsom ngerm nnd ntot] = count_val_results(val)

nsom = sum(strncmpi(val,'som',3));
ngerm = sum(strncmpi(val,'germ',4));
nnd = sum(strncmpi(val,'not',3));
ntot = nsom+ngerm+nnd;
