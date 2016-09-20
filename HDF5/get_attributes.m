function attr = get_attributes(D,field)
%GET_ATTRIBUTES get the attributes associated with FIELD in datastruct D.
%
%

%---
% $Id$
% $Date: 2007-11-30 13:53:24 -0500 (Fri, 30 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$


idx = strmatch(field,D.fieldnames,'exact');

attr = D.fielddata{idx}.Attributes;
