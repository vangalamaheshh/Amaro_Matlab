function indices  = findsupdat(D,supacc,equalto,relate)
% FINDSUPDAT return indices of D matching specified value.
%
%   INDICES = FINDSUPDAT(D,SUPACC,EQUALTO) gets the column indices of D for
%   which D.supdat(SUPACC,:) == EQUALTO.  If EQUALTO is not defined, Gets
%   the column indices of D for which SUPACC ~= 0. D.supdat may be 2 or 3
%   dimensional.  If D.supdat is 3 dimensional, EQUALTO is compared to the
%   dim 3 sum of supdat.  Optional input RELATE sets the way 3 dimensional
%   supdat will be collaped into 2 dimensions for comparison.  Default is
%   'and' (each .supdat value must be true for a given sample if that
%   samples index is to be returned.)
%
%       Revisions:
%           12 Oct 07 -- Created by Jen Dobson (jdobson@broad.mit.edu)
%
%---
% $Id$
% $Date: 2007-11-19 17:28:20 -0500 (Mon, 19 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$

if ndims(D.supdat) == 2
    
    if exist('equalto','var')
        indices = find(D.supdat(find_supid(D,supacc),:)==equalto);
    else
        indices = find(D.supdat(find_supid(D,supacc),:)~=0);
    end
    
else
    
    if ~exist('relate','var')
        relate = 'and'
    end
    
    switch relate
        case 'sum'
            
           dat = sum(D.supdat(find_supid(D,supacc),:,:),3);
           
        case 'mean'
            
            dat = mean(D.supdat(find_supid(D,supacc),:,:),3);
           
        case 'or'
            
            dat = zeros(1,size(D.supdat,2));
            for k = 1:size(D.supdat,3)
                dat = equalto.*or(dat, D.supdat(find_supid(D,supacc),:,k)==equalto);
            end
            
        case 'and'
            
            dat = ones(1,size(D.supdat,2));
            for k = 1:size(D.supdat,3)
                dat = equalto.*and(dat, D.supdat(find_supid(D,supacc),:,k)==equalto);
            end

    end
    
    if exist('equalto','var')
        indices = find(dat==equalto);
    else
        indices = find(dat~=0);
    end
    
end

    