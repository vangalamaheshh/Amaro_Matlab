function print_pdf(fname)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

cf=num2str(gcf);
rend=get(gcf,'renderer');
set(gcf,'renderer','none');
print(['-f' cf],['-dpsc2'],[ fname '.ps' ]);
set(gcf,'renderer',rend);

ps2pdf([ fname '.ps'],[fname '.pdf']);

