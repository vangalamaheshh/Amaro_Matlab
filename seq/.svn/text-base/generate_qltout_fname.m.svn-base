function fname=generate_qltout_fname(i,freeze)
if ~exist('freeze','var') || isempty(freeze), freeze=default_freeze; end

lanes = load_lanelist;
freezes = load_freezelist;

if freeze<1 || freeze>slength(freezes), error('No such freeze'); end

idx = find(lanes.FLE==i);
if length(idx)~=1
  fprintf('Lane index not found in lanelist!\n');
  fname = 'none';
else
  if strcmp(lanes.tumornormal{idx},'tumor'), fname = freezes.tumorbasedir{freeze};
  else fname = freezes.normalbasedir{freeze}; end
  fname = [fname '/alignments'];
%%  fname = ['/home/radon01/kiran/tmp/tcga/tcga-' lanes.tumornormal{idx} '-freeze' ...
%%    num2str(freeze) '-' freezes.date{freeze} '/alignments'];

%  if strcmp(lanes.TN{idx},'N')
%    fname = '/home/radon01/kiran/tmp/tcga/tcga-normal-freeze1-081121/alignments';
%  else
%    fname = '/home/radon01/kiran/tmp/tcga/tcga-tumor-freeze1-081121/alignments';
%  end

  fname = [fname '/' lanes.FC{idx} '.' lanes.lane{idx}];

  if lanes.has_pair(idx)
    fname = [fname '/pairfinder/' lanes.FC{idx} '.' lanes.lane{idx} '.end' lanes.end{idx} '.qltout'];
  else 
    if strcmp(lanes.date{idx},'None')
       fprintf('This lane was single-endread-only, and it lacks that endread.\n');
       fname = 'none';
    else
       fname = [fname '/' lanes.date{idx} '/' lanes.FC{idx} '.' lanes.lane{idx} '.qltout'];
    end
  end
end
