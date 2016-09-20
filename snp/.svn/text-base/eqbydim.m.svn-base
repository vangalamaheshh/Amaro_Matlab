function eqarray = eqbydim(A,dim)
%EQBYDIM checks for equality of array A along dimension DIM.
%
%       EQARRAY = EQBYDIM(A,DIM) returns N-1 dimensional array EQARRAY. An
%       element of EQARRAY is 1 if all elements in the vector along
%       dimension DIM are equal.  Otherwise the element is 0.
%       
%       Example:  A = [1 4 7 2 3; 1 3 7 4 3; 1 4 7 5 3]
%                 E = eqbydim(A,1);
%                 E = [1 0 1 0 1];


if dim > ndims(A)
    error('Dimension of input matrix A is greater than input dimension dim');
end

sizeA = size(A);
savedims = setdiff([1:ndims(A)],dim);



Aro = permute(A,[dim savedims]);

eqarray = ones(size(Aro(1,:)));
for k = 2:size(Aro,1)
    eqarray = eqarray & (Aro(k-1,:) == Aro(k,:));
end


eqarray = reshape(eqarray,sizeA(savedims));

