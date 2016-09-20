function [supacc, supdesc, supdat, mergsupdat] = preproc_mergesup(Cin,mkuarray)
%PREPROC_MERGESUP merge .sup* when merging data structures in Cin.
%
%   [SUPACC,SUPDESC,SUPDAT] = PREPROC_MERGESUP(CIN,sortarray) is a
%   subroutine to PREPROC_MERGEPLATFORMS.  CIN is a cell array of data
%   structures that are being merged, SORTIDX is the M X N reordering index
%   array for the samples in Cin.  (M is number of samples, N is length of
%   CIN).  The new SUPACC and SUPDESC are just the .supacc and .supdesc of
%   the first data struct in Cin  (Cin{1}).  SUPDAT is the new supdat with
%   dims the same as Cin{1}.supdat.  MERGSUPDAT is a 3 dimensional array
%   where each increment in dimension 3 gives the supdat for one of the
%   platforms.   (i.e. MERGSUPDAT is a T X S X P array where T is number of locations, S is
%   number of samples, P is number of platforms.  
%
%   The merge rules for SUPDAT are: 
%       *array specific information merges to NaNs; 
%
%       *include_ fields (listing the samples to include for the various
%       steps of pipeline) merge to the maximum number across dimension 3;
%
%  NEED TO ADD: what if both platforms don't have the same supaccs??
%

%           Revisions:
%               -- 17 Oct 07:  Function created by Jen Dobson
%               (jdobson@broad.mit.edu).
%
%               -- 1 Nov 07:  supdat will be one - dimensional.  Fields
%               that are array specific merge to NaNs.  mergedsupdat will
%               give multidimensional original supdat info.
%
%               -- 7 Nov 07:  Updated to check for error
%       
%---
% $Id$
% $Date: 2007-11-19 17:28:20 -0500 (Mon, 19 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$



supacc = Cin{1}.supacc;

supdesc = Cin{1}.supdesc;

tmpdat = zeros(size(Cin{1}.supdat,1),size(mkuarray{1},1),length(Cin));

for k = 1:length(Cin)
    [matx, idx1, idx2 ] = match_string_sets(supacc,Cin{k}.supacc);
    tmpdat(:,:,k) = Cin{k}.supdat(idx1,mkuarray{k});
end

mergsupdat = tmpdat;



arrayspecific = {'array','dup','core','good','loh','batch','rep','force'} ;
[m,dum,asidx] =  match_string_sets(arrayspecific,lower(strtok(cellstr(supacc),':')));
tmpdat(asidx,:,:) = Inf;   %fill the array specific info with Inf  (since when we check for matches Inf == Inf)

%% Check for mismatches in array list include dat;  make mismatched data take on highest value (?)


% Check include_* supdat
aliidx = strmatch('include_',supacc);  %find the supdat rows corresponding to include_*

if ~isempty(aliidx)
    [badali,badalisamp] = find(~eqbydim(tmpdat(aliidx,:,:),3));  %find the samples that don't have the same include_* label in all merged samples

    if ~isempty(badali)
        sampnames = {Cin{1}.sis(mkuarray{1}).uniqid};
        for kk = 1:length(badali)
            warning('Mismatched array list include data for sample: %s in column: %s \n  Replacing with max values.',...
                char(sampnames(badalisamp(kk))),supacc(aliidx(badali(kk)),:));

        end
    end
end

for z = 1:size(tmpdat,3)
    tmpdat(aliidx,:,z) = max(tmpdat(aliidx,:,:),[],3);
end


matched = ones(size(tmpdat,1),size(tmpdat,2));

for i = 2:size(tmpdat,3)
    matched = matched & (tmpdat(:,:,i) == tmpdat(:,:,i-1));
end

supdat = zeros(size(tmpdat,1),size(tmpdat,2));
supdat = tmpdat(:,:,1) .* matched;
supdat(find(~matched)) = NaN;

[nmi,nmj] = find(isnan(supdat));

unmi=unique(nmi);

if ~isempty(unmi)
    for ii = unmi'
        badsamps = find(nmi == ii)';
        warning(['Mismatched supdat for %s in columns:  ' num2str(badsamps) '\nReplacing with NaNs.'],deblank(char(supacc(ii,:))));
    end
end

supdat(isinf(supdat)) = NaN;

    
