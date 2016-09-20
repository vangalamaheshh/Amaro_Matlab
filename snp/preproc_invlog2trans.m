function Cnew = preproc_invlog2trans(Craw,field)
%PREPROC_INVLOG2TRANS take inverse log of raw data.
%
%CNEW = PREPROC_INVLOG2TRANS(CRAW) returns the data structure CNEW 
% after inverse log2 transform.
%
% Optional FIELD input allows selection of different field to transform.
%
%       History:
%           26 Sept 07:  Created by Jen Dobson (jdobson@broad.mit.edu)
%           10 Jan 08: Added inner loops for memory conservation.
%           (jdobson@broad.mit.edu)
%---
%  $Id$
%  $Date: 2008-02-01 17:31:09 -0500 (Fri, 01 Feb 2008) $
%  $LastChangedBy: jdobson $
%  $Rev$



if ~iscell(Craw)
    Ctemp = {Craw};
    Craw = Ctemp;
    clear Ctemp;
    wascell = 0;
else
    wascell = 1;
end


    if ~exist('field','var')
     field = 'dat';
end


for k = 1:length(Craw)

    chunkdims = getmemchunkdims(Craw{k},field);

    datsize = getsize(Craw{k},field);

    ii = 0;
    jj = 0;

    while ii < datsize(1)
        thisloopdim1idx = (ii+1):min(ii+chunkdims(1),datsize(1));
        while jj < datsize(2)
            thisloopdim2idx = (jj+1):min(jj+chunkdims(2),datsize(2));
            Craw{k}.(field)(thisloopdim1idx,thisloopdim2idx) = 2.^(Craw{k}.(field)(thisloopdim1idx,thisloopdim2idx)+1);
            jj = jj + chunkdims(2);
        end
        ii = ii + chunkdims(1);
        jj = 0;
    end
    Craw{k} = add_history(Craw{k},mfilename);
end

Cnew = Craw;
if ~wascell
    Cnew = Cnew{1};
end
