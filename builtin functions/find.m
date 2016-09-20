function [varargout] = find(varargin)
%FIND   Find indices of nonzero elements.
%   I = FIND(X) returns the linear indices of the array X that are nonzero.
%   X may be a logical expression. Use IND2SUB(I,SIZE(X)) to calculate
%   multiple subscripts from the linear indices I.
% 
%   I = FIND(X,K) returns at most the first K indices of X that are
%   nonzero. K must be a positive integer, but can be of any numeric type.
%
%   I = FIND(X,K,'first') is the same as I = FIND(X,K).
%
%   I = FIND(X,K,'last') returns at most the last K indices of X that are
%   nonzero.
%
%   [I,J] = FIND(X,...) returns the row and column indices instead of
%   linear indices into X. This syntax is especially useful when working
%   with sparse matrices.  If X is an N-dimensional array where N > 2, then
%   J is a linear index over the N-1 trailing dimensions of X.
%
%   [I,J,V] = FIND(X,...) also returns a vector V containing the values
%   that correspond to the row and column indices I and J.  If X is a
%   logical expression, then V will contain the values returned after
%   evaluating that expression.
%
%   Example:
%      A = magic(3)
%      find(A > 5)
%
%   finds the linear indices of the 4 entries of the matrix A that are
%   greater than 5.
%
%      [rows,cols,vals] = find(speye(5))
%
%   finds the row and column indices and nonzero values of the 5-by-5
%   sparse identity matrix.
%
%   See also SPARSE, IND2SUB, RELOP, NONZEROS.

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 5.11.4.7 $  $Date: 2004/12/06 16:34:08 $
%   Built-in function.

if nargout == 0
  builtin('find', varargin{:});
else
  [varargout{1:nargout}] = builtin('find', varargin{:});
end
