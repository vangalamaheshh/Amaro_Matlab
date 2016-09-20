function set_verbose_file(filename)
% SET_VERBOSE_FILE(FILENAME) sets the global VERBOSE_FILE  to FILE.
%
%    If set, all verbose output will be written to FILENAME
%
%    See also VERBOSE, SET_VERBOSE_LEVEL

%---
% $Id$
% $Date: 2008-01-10 11:26:15 -0500 (Thu, 10 Jan 2008) $
% $LastChangedBy: jdobson $
% $Rev$ 



if ~is_full_filename(filename)
    filename = [pwd filesep filename];
    warning(['Changing filename to full file path: ' filename]);
end


global VERBOSE_FILE

VERBOSE_FILE=filename;
