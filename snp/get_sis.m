function sisdata = get_sis(D,sisfield,sisidx)
%GET_SIS get the contents of a field of the D's sis returned in a cell
%array.
%
%       SISDATA = GET_SIS(D,SISFIELD,SISIDX)
%
%   Note: Written to abstract away from D implementation as struct or
%   datastruct:
%            can't make datastruct object return a comma-separated
%                  list (bug in matlab; no work-around.  see:
%                  http://www.mathworks.com/support/solutions/data/1-1AB
%                  OD.html?solution=1-1ABOD
%
%       Revisions:
%           12 Dec 07: function created (jdobson)

%---
% $Id$
% $Date: 2008-02-14 09:08:27 -0500 (Thu, 14 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$


if ~isfield(D,'sis')
    error('No sis field found')
end


if isa(D,'struct')

    if ~exist('sisidx','var')

        sisdata = {D.sis.(sisfield)};

    else
        sisdata = {D.sis(sisidx).(sisfield)};
    end


elseif isa(D,'datastruct')

    if ~exist('sisidx','var')

        sisdata = D.sis.(sisfield);

    else

        theSIS = D.sis;
      
        sisdata = {theSIS(sisidx).(sisfield)};
    end

else

    error('GET_SIS first input argument must be of type STRUCT or DATASTRUCT');

end

