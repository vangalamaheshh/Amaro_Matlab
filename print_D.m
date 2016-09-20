function print_D(fname,formats,landscape,papersize,use_renderer,matlab_version)
%
% ---
% $Id$
% $Date: 2007-09-17 14:46:20 -0400 (Mon, 17 Sep 2007) $
% $LastChangedBy: rameen $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

 
  
if ~exist('landscape','var')
  landscape=0;
end

if ~exist('matlab_version','var') || isempty(matlab_version)
  matlab_version = '-v7.3';
end

if landscape
  set(gcf,'PaperOrientation','landscape');
  pp=get(gcf,'PaperPosition');
  set(gcf,'PaperPosition',[ 0.25 0.25 10.25 pp(4)*10.25/pp(3)]); 
end
rend=get(gcf,'renderer');
if ~exist('use_renderer','var') || isempty(use_renderer) || ~use_renderer 
  set(gcf,'renderer','none');
end
cf=num2str(gcf);
for i=1:length(formats)
  cur=formats{i};
  if ~isempty(find(fname=='.')) && ~strcmp(cur{1},'pdf')
    fname_i=[fname '.' cur{1}];
  else
    fname_i=fname;
  end
%  disp(cur{1});
  if strcmp(cur{1},'pdf')
    print_pdf(fname_i);
  elseif strcmp(cur{1},'fig')
    hgsave([fname_i],matlab_version)
  else
    if length(cur)==1
      print(['-f' cf],['-d' cur{1}],[ fname_i ]);
      % add a piece to the name which describes the parameters used
      % (reduced,...)
    else
      print(['-f' cf],['-d' cur{1}],cur{2:end},[ fname_i ]); %was [fname '.' cur{1}]
    end  
  end
end
rend='painters';
set(gcf,'renderer',rend);
disp(['setting renderer to ' rend]);

%print('-f1','-dpng','-r180',[ fname add_dist_st '.png'])
%print('-f1','-depsc',[ fname add_dist_st '.eps']);
%print('-f1','-dpdf',[ fname add_dist_st '.pdf']);




