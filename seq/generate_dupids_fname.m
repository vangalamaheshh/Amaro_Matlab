function fname=generate_dupids_fname(i,freeze)
if ~exist('freeze','var') || isempty(freeze), freeze=2; end

lanes = load_lanelist;
freezes = load_freezelist;

if freeze<1 || freeze>slength(freezes), error('No such freeze'); end

idx = find(lanes.FLE==i);
if length(idx)~=1
  fprintf('Lane index not found in lanelist!');
  fname = 'none';
else
  if ~lanes.has_pair(idx)
    fprintf('This lane was single-endread-only, and such could not be de-duped.');
    fname = 'none';
  else
    if strcmp(lanes.tumornormal{idx},'tumor'), fname = freezes.tumorbasedir{freeze};
    else fname = freezes.normalbasedir{freeze}; end

%%  fname = ['/home/radon01/kiran/tmp/tcga/tcga-' lanes.tumornormal{idx} '-freeze' ...
%%    num2str(freeze) '-' freezes.date{freeze} '/'];

%    if strcmp(lanes.TN{idx},'N')
%      fname = '/home/radon01/kiran/tmp/tcga/tcga-normal-freeze1-081121/';
%    else
%      fname = '/home/radon01/kiran/tmp/tcga/tcga-tumor-freeze1-081121/';
%    end
    fname = [fname '/meta/wg.e4/' lanes.FC{idx} '.' lanes.lane{idx}...
       '.end' lanes.end{idx} '.qltout.dupids'];
  end
end
