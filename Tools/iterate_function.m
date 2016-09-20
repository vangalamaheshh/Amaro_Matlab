function D = iterate_function(D,field,fulldim,fh,varargin)
% D = ITERATE_FUNCTION(D,FIELD,FULLDIM,FUNCTION_HANDLE,VARARGIN)

if ~exist('fulldim','var')
    fulldim = [];
end


chunkdims = getmemchunkdims(D,field,fulldim);

datsize = getsize(D,field);

ii = 0;
jj = 0;



while ii < datsize(1)
        thisloopdim1idx = (ii+1):min(ii+chunkdims(1),datsize(1));
        while jj < datsize(2)
            thisloopdim2idx = (jj+1):min(jj+chunkdims(2),datsize(2));
            D.(field)(thisloopdim1idx,thisloopdim2idx) = fh(D.(char)(thisloopdim1idx,thisloopdim2idx),varargin);
            jj = jj + chunkdims(2);
            fprintf(1,'.');
        end
        ii = ii + chunkdims(1);
        jj = 0;
    end
 
end
