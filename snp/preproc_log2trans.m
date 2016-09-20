function Cnew = preproc_log2trans(Craw,replacesmallvals,floorval,field)
%PREPROC_LOG2TRANS take log of raw data.
%
%Cnew = PREPROC_LOG2TRANS(CRAW,REPLACESMALLVALS,FLOORVAL,field) returns the data
%structure CRAW in CNEW after log2 transformation OF CRAW.dat.  An optional
%toggle set by REPLACESMALLVALS=1 will replace values of .dat less than
%FLOORVAL with FLOORVAL.  CRAW can be a gistic data structure or a cell
%array of gistic data structures.  The default value of FLOORVAL is 1.
%
%  Optional field input allows you to select a different field to transform.
%
%
%       History:
%           26 Sept 07:  Created by Jen Dobson (jdobson@broad.mit.edu)
%           10 Jan 08:  Added inner loops for memory conservation.
%           (jdobson@broad.mit.edu).
%
%---
%  $Id$
%  $Date: 2008-06-19 15:10:18 -0400 (Thu, 19 Jun 2008) $
%  $LastChangedBy: jdobson $
%  $Rev$

if ~exist('replacesmallvals','var')
    replacesmallvals = 1;
end


if ~exist('floorval','var')
    floorval = 1;
end


if replacesmallvals
    verbose(['Replacing values <' num2str(floorval) ' with ' num2str(floorval)],30);
end




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

    if replacesmallvals
        repidx = find(Craw{k}.(field) < floorval);
        Craw{k}.(field)(repidx)=floorval;  %replace small values with floor value
    end
    
    
      chunkdims = getmemchunkdims(Craw{k},field);

    datsize = getsize(Craw{k},field);

    ii = 0;
    jj = 0;
    
    verbose(['Taking log2 of intensities by chunks: ' num2str(chunkdims)],30)

    while ii < datsize(1)
        thisloopdim1idx = (ii+1):min(ii+chunkdims(1),datsize(1));
        while jj < datsize(2)
            thisloopdim2idx = (jj+1):min(jj+chunkdims(2),datsize(2));
            Craw{k}.(field)(thisloopdim1idx,thisloopdim2idx) = log2(Craw{k}.(field)(thisloopdim1idx,thisloopdim2idx))-1;
            jj = jj + chunkdims(2);
        end
        ii = ii + chunkdims(1);
        jj = 0;
    end
 
    verbose(['Log2 complete'],30);
    Craw{k} = add_history(Craw{k},mfilename,floorval);
end



Cnew = Craw;
if ~wascell
    Cnew = Cnew{1};
end
