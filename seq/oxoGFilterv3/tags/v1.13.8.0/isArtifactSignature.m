function result = isArtifactSignature(refAlleles, altAlleles, altSignature)
%
% result = isArtifactSignature(refAlleles, altAlleles)
% 
% refAlleles and altAlleles must be vectors of single characters in {'C', 'G', 'A', 'T', 'N'}
%   These should be in cells and ref and alt should be of the same length
%   and corresponding.
%
%   lower case letters are acceptable
%
%   returns 0 or 1 (false and true, respectively)

if length(refAlleles) ~= length(altAlleles)
   error('ref and alt vectors are not the same length') 
end

if nargin < 3 || isempty(altSignature)
   altSignature = 'cagt'; 
end

if length(altSignature) ~= 4
   error('Currently only two mutations are supported at one time.  One can be supported with duplication (e.g. caca -- will only look at C>A mutations)') 
end

r = lower(refAlleles);
a = lower(altAlleles);

result = (strcmp(r,altSignature(1)) & strcmp(a, altSignature(2)));
result = result | (strcmp(r,altSignature(3)) & strcmp(a, altSignature(4)));

% Check for DNP/TNPs
refAlleleLengths = cellfun(@length,r);
altAlleleLengths = cellfun(@length,a);

if ~all(refAlleleLengths == altAlleleLengths)
   error('Not all ref and alt allele lengths are the same.  Is this file exclusively SNP, DNP, TNP?') 
end

if any(refAlleleLengths > 1)
   indices = find(refAlleleLengths > 1);
   warning(['DNP or TNP detected.  These have not been evaluated.  Indices: ' num2str(indices')]) 
end

