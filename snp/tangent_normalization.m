function [Y,Npi,weights] = tangent_normalization(X,N,Npi,varargin)
%TANGENT_NORMALIZATION low level tangent normalization routine
%
%   [Y,NPI] = tangent_normalization(X,N,NPI)
%
% Input X is a column vector to normalize and input N is an array
% of normal vectors. Returns Y, a normalized vector formed by subtracting
% the projection of X onto the plane defined by linear combinations of the
% vectors of N. The input/output variable NPI is the pseudoinverse of N,
% which can be kept in memory to avoid (time-consuming) recalculation.

% unless provided by input, generate the pseudoinverse
if ~exist('Npi','var') || isempty(Npi)
    Npi = pinv(N);
end
% subtract input's projection on normal plane from input
weights = Npi*X;
Y = X - N*weights;
