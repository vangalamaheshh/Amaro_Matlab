function output = itrfcn1(D,field,fulldim,fh,varargin)
% D = ITERATE_FUNCTION(D,FIELD,FULLDIM,FUNCTION_HANDLE,VARARGIN)

if ~exist('fulldim','var')
    fulldim = [];
end


chunkdims = getmemchunkdims(D,field,fulldim);

datsize = getsize(D,field);

ii = 0;
jj = 0;


iterdim = find(chunkdims ~= datsize);


if ~isempty(iterdim)

output = zeros(1,datsize(iterdim));
    if iterdim == 1
        output = output';
    end

    fprintf(1,'Iterating function:');
    while ii < datsize(iterdim)
        thisloopidx = (ii+1):min(ii+chunkdims(iterdim),datsize(iterdim));
        if iterdim == 1
            output(thisloopidx) = fh(D.(field)(thisloopidx,:),varargin{:});
        elseif iterdim == 2
            output(thisloopidx) = fh(D.(field)(:,thisloopidx),varargin{:});
        end

        fprintf(1,'.');
        ii = ii + chunkdims(iterdim);
      
       
    end

    fprintf(1,'\n');

else
    output = fh(D.(field),varargin{:});
end
