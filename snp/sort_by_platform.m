function M = sort_by_platform(D)
%
%DEAL_BY_PLATFORM Sort cells of data structure by platform and deal to new
%data structure.
%
%   M = DEAL_BY_PLATFORM(D) assigns the cell contents of D to M, such that 
%   each cell of M represents a unique platform.  The data structures in 
%   the cells of D must contain platform information in the .sis field, and 
%   each cell of D should have data from only one platform.  The standard
%   platform identifiers in D.sis.platform are strings '100X', '100H',...
%
%   History:
%           20 Sept 07 -- Created by Jen Dobson (jdobson@broad.mit.edu)
%


fun = @(X) getfield(getfield(X,'sis'),'platform');  %function handle for cellfun
platforms = cellfun(fun,D,'UniformOutput',0);  %get the platform from each cell

[uplats, PI, UPJ] = unique(platforms);

gun = @(X) find(X==UPJ);
DidxM = cellfun(gun,num2cell([1:length(uplats)]),'UniformOutput',0);

hun = @(X) D(X);
M = cellfun(hun,DidxM,'UniformOutput',0);

