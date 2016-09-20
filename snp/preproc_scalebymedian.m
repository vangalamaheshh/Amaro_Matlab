function [Cnew] = preproc_scalebymedian(Craw,dataislog)
%PREPROC_SCALEBYMEDIAN Scale data by median value
%
%   [CNEW] = PREPROC_SCALEBYMEDIAN(CRAW,DATAISLOG) scales each column
%   of sample data in struct CRAW.dat to the data column's median value and
%   returns the scaled data in CNEW. DATAISLOG is used to specify whether
%   the input data is log transformed (1) or linear (0). Default for
%   DATAISLOG is 1.  CRAW can be a struct or a cell array of structures.  
%---
% $Id$
% $Date: 2008-01-17 10:33:36 -0500 (Thu, 17 Jan 2008) $
% $LastChangedBy: jdobson $
% $Rev$


    


if ~exist('dataislog','var')
    dataislog = 1;
end

if ~iscell(Craw)
    Cnew = {Craw};
    wascell = 0;
else
    Cnew = Craw;
    wascell = 1;
end

verbose('Normalizing samples to median values',20);

for k = 1:length(Cnew)

    memchunks = getmemchunkdims(Cnew{k},'dat',1);
    maxsamps = memchunks(2);
verbose(['Chunk size: ' num2str(memchunks)],30);
    jj = 0;

    while jj < size(Cnew{k}.dat,2)
        thisendidx = min(jj+maxsamps,size(Cnew{k}.dat,2));
        dat = Cnew{k}.dat(:,jj+1:thisendidx);
        med = median(dat,1);
        if dataislog
            dat = dat-repmat(med,size(dat,1),1);
        else
            dat = dat./repmat(med,size(dat,1),1);
        end
        Cnew{k}.dat(:,jj+1:thisendidx)=dat;
        jj = jj+maxsamps;
    end

end

if ~wascell
    Cnew = Cnew{1};
end


