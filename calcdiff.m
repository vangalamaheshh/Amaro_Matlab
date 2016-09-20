function D = calcdiff(D,field,N,dim,difffield)
% D = CALCDIFF(D,FIELD,N,DIM,DIFFFIELD)
%
%  D is a datastructure, FIELD is the field to find the difference of, N
%  is the order of the diff calculation, DIM is the dimension along which
%  to take the diff, DIFFFIELD is optional and is the fieldname of the
%  difference matrix.  (jdobson@broad.mit.edu)
%

%  ---
% $Id$
% $Date: 2008-02-28 15:30:40 -0500 (Thu, 28 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$

if ~exist('difffield','var') || isempty(difffield)
    difffield = 'diff';
end


if strcmp('datastruct',class(D))


    newfname = strcat(regexp(get_datafile(D,field),[regexp_filesep '.+' regexp_filesep],'match'),difffield,'.h5');
    newfname = get_new_filenames(newfname);
    dims = getsize(D,field);
    deltadim = [0 0];
    deltadim(dim) = 1;
    dims = dims - deltadim;
    D = add_diskfield(D,newfname,difffield,dims,'single');


    chunkdims = getmemchunkdims(D,difffield,1);


    ii = 0;
    jj = 0;

    if ~isequal(chunkdims,dims)


        fprintf(1,'Iterating diff:');


        while ii < dims(1)
            thisloopdim1idx = (ii+1):min(ii+chunkdims(1),dims(1));
            while jj < dims(2)
                thisloopdim2idx = (jj+1):min(jj+chunkdims(2),dims(2));

                D.(difffield)(thisloopdim1idx,thisloopdim2idx) = diff(D.(field)...
                    (union(thisloopdim1idx,thisloopdim1idx(end) + deltadim(1)),union(thisloopdim2idx,thisloopdim2idx(end) + deltadim(2))),...
                    N,dim);


                fprintf(1,'.');
                jj = jj + chunkdims(2);
            end

            ii = ii + chunkdims(1);
            jj = 0;

        end

        fprintf(1,'\n');

    else
        D.(difffield) = diff(D.(field),N,dim);

    end

else
    D.(difffield) = diff(D.(field),N,dim);
end

