function fname=generate_illogicals_fname(i,freeze)
if ~exist('freeze','var') || isempty(freeze), freeze=default_freeze; end

lanes = load_lanelist;
freezes = load_freezelist;

if freeze~=-1 && (freeze<1 || freeze>slength(freezes)), error('No such freeze'); end

idx = find(lanes.FLE==i);
if length(idx)~=1
  fprintf('Lane index not found in lanelist!');
  fname = 'none';
else
  if freeze==-1
    fname = '/wga/scr1/cancer_Illumina_WGS/illogicals_tmp';
  else
    if strcmp(lanes.tumornormal{idx},'tumor'), fname = freezes.tumorbasedir{freeze};
    else fname = freezes.normalbasedir{freeze}; end
    fname = [fname '/alignments/' lanes.FC{idx} '.' lanes.lane{idx} '/pairfinder'];
%%    fname = ['/home/radon01/kiran/tmp/tcga/tcga-' lanes.tumornormal{idx} '-freeze' ...
%%    num2str(freeze) '-' freezes.date{freeze} '/alignments/' ...
%%    lanes.FC{idx} '.' lanes.lane{idx} '/pairfinder'];
  end

%  if strcmp(lanes.TN{idx},'N')
%    fname = '/home/radon01/kiran/tmp/tcga/tcga-normal-freeze1-081121/alignments';
%  else
%    fname = '/home/radon01/kiran/tmp/tcga/tcga-tumor-freeze1-081121/alignments';
%  end

  if lanes.has_pair(idx)
    fname = [fname '/' lanes.FC{idx} '.' lanes.lane{idx} '.end' lanes.end{idx} '.qltout.untrusted.illogicals'];
  else 
    fprintf('This lane was single-endread-only, and so it lacks an illogicals file.');
    fname = 'none';
  end
end
