function D = itrfcn2(D,field,newfield,fulldim,fh,varargin)
% D = ITRFCN2(D,FIELD,NEWFIELD,FULLDIM,FUNCTION_HANDLE,VARARGIN)

if ~exist('fulldim','var')
    fulldim = [];
end


chunkdims = getmemchunkdims(D,field,fulldim);

datsize = getsize(D,field);

ii = 0;
jj = 0;



if ~isequal(chunkdims,datsize)
    if ~isfield(D,newfield)
    D.(newfield) = zeros(datsize(1),datsize(2));
    end
    
    fprintf(1,'Iterating function:');
    
    
    while ii < datsize(1)
        thisloopdim1idx = (ii+1):min(ii+chunkdims(1),datsize(1));
        while jj < datsize(2)
            thisloopdim2idx = (jj+1):min(jj+chunkdims(2),datsize(2));
        
            D.(newfield)(thisloopdim1idx,thisloopdim2idx) = fh(D.(field)(thisloopdim1idx,thisloopdim2idx),varargin{:});
      

        fprintf(1,'.');
        jj = jj + chunkdims(2);
        end
        
        ii = ii + chunkdims(1);
        jj = 0;
       
    end
    
    fprintf(1,'\n');

else
    D.(newfield) = fh(D.(field),varargin{:});
end


