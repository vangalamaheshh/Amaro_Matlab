function chunkdims = getmemchunkdims(D,field,fulldim)
%GETMEMCHUNKDIMS get the dimensions of the maximal part of D.FIELD
%allowed in memory as set by D's MemLimit.
%
%   CHUNKDIMS = GETCHUNKDIMS(D,FIELD,FULLDIM)
%           FULLDIM: Default is to take full vector along the shorter
%           dimension.  If FULLDIM is specified, take a full vector along
%           dimension FULLDIM. 

%       Revisions:
%           * Function added: 10 Jan 07
%
%---
% $Id$
% $Date: 2008-02-28 17:08:14 -0500 (Thu, 28 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$



%% Get size of one element of the array to figure out how big
%the read chunks can be

if ~isdiskfield(D,field)
    chunkdims = getsize(D,field);
else
    mlimit = getmetadata(D,'MemLimit');
    idx = strmatch(field,D.fieldnames,'exact');
    dims = getsize(D,field);

    if isempty(mlimit)
        chunkdims = dims;
    else
        sampleelement = cast(0,htype_to_mtype(D.fielddata{idx}.Datatype.Class)); %#ok
        wh = whos('sampleelement');
        eltsize = wh.bytes;






        %% Set the chunksize by first taking as much as possible along the smaller
        % dimension (MINDIMIDX) (hopefully the whole thing); then take as much as
        % possible along the larger dimension (MAXDIMIDX)
        chunkdims = zeros(1,2);

        if ~exist('fulldim','var') || isempty(fulldim) || ~(fulldim == 1 || fulldim == 2)
            [fulldimsize,fulldimidx] = max(dims);
            [partdimsize,partdimidx] = min(dims);
        else
            fulldimsize = dims(fulldim);   %if fulldim is specified, treat it as the smaller dimension
            fulldimidx = fulldim;
            partdimidx = setdiff([1 2],fulldim);
            partdimsize = dims(partdimidx);
        end

        vectsize = eltsize.*fulldimsize;
        if vectsize > mlimit
            chunkdims(fulldimidx) = floor(mlimit./eltsize);
            chunkdims(partdimidx) = 1;
        else
            chunkdims(fulldimidx) = fulldimsize;
            chunkdims(partdimidx) = min(floor(mlimit./vectsize),partdimsize);

        end

    end
end