function ps2pdf(psname,pdfname,dont_remove_ps)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('dont_remove_ps','var')
  dont_remove_ps=0;
end

if strcmp(computer,'GLNXA64')
  st='/broad/tools/apps/matlab73/';
%  st='/util/matlab7.2/';
else
%  st='/ibm_local/util/matlab7.2/';
  st='/broad/tools/apps/matlab73/';
end

unix([st 'sys/ghostscript/bin/' lower(computer) '/gs' ...
      ' -dNOPAUSE' ...
      ' -I"' st 'sys/ghostscript/ps_files"' ...
      ' -I"' st 'sys/ghostscript/fonts"' ...
      ' -sDEVICE=pdfwrite -sOutputFile=' pdfname ' -dBATCH ' psname ]);
if ~dont_remove_ps && exist(pdfname,'file')
  unix(['rm ' psname]); 
end
