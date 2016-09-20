function lego(z,P)
% input should be a 8x12 matrix of heights for the "LEGO" plot
% tries to convert other formats

if ~exist('P','var'), P=[]; end

if iscell(z), legos(z,P); return; end
if isstruct(z) && isfield(z,'context65'), legomaf(z,P); return; end
if length(size(z))==2 & all(size(z)==[1 96]), z = factor2lego(z); end
if length(size(z))==2 & all(size(z)==[96 1]), z = factor2lego(z'); end
if length(size(z))==4 & all(size(z)==[4 4 4 4] | size(z)==[4 2 4 4]), z = matrix2lego(z); end
if ~(length(size(z))==2 & all(size(z)==[8 12])), error('unknown format'); end
 
P = impose_default_value(P,'imagesc',false);  % if true, displays as imagesc instead of bar3_with_colors
P = impose_default_value(P,'zmax',[]);
if P.imagesc
  if ~isempty(P.zmax), z = min(z,P.zmax); end
  imagesc(z);
  pp = {'linewidth',2,'color',[0 0 0]};
  line([4.5 4.5],ylim,pp{:});line([8.5 8.5],ylim,pp{:});line(xlim,[4.5 4.5],pp{:});
else
  c = get_LEGO_colors();
%  z3 = repmat(z,[1 1 3]); c(isnan(z3))=1; % white   %%% for some reason this is causing bizarre color scrambling in context of legos.m
  bar3_with_colors(z,c);
  if ~isempty(P.zmax), zlim([0 P.zmax]); end
end
set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);



