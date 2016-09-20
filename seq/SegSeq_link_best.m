function SegSeq_link_best(sample,P)
% For a sample that SegSeq has been run on,
%   goes into the output directory and checks each .seg.txt file,
%   sees how many segments were generated.
% Flowcell with the fewest segments is judged "best".
% Links are created in the "lawrence" directory structure,
%   "SegSeq_results.txt" and "SegSeq_results.png",
%   pointing to the output from the best flowcell.
%
% Mike Lawrence 2009-09-23

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end
P = impose_default_value(P,'force_use_flowcell',[]);

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Multiple samples not supported'); end
if ~contains(sample,'wgs') || contains(sample,'-') || ~contains(sample,'/')
  error('<sample> should be of form ov/0751/wgs');
end
name2 = upper(regexprep(sample,'/','-'));
name2 = regexprep(name2,'-WGS$','');

lawdir = ['/xchip/tcga_scratch/lawrence/' sample];
ssdir = ['/xchip/tcga_scratch2/ng/' name2 '/matfiles'];

if ~exist(ssdir,'dir'), error('Not found: %s',ssdir); end
ssd = dir([ssdir '/*.seg.txt']);
if length(ssd)==0, error('No seg files in %s',ssdir); end

F = [];
for i=1:length(ssd)
  F.segname{i,1} = [ssdir '/' ssd(i).name];
  x = load_lines(F.segname{i});
  F.nseg(i,1) = length(x)-1;
  F.pngname{i,1} = regexprep(F.segname{i},'.seg.txt$','.png');
  d = dir(F.pngname{i});
  F.haspng(i,1) = ~isempty(d);
end

if ~any(F.haspng), error('No png files in %s',ssdir); end
F = reorder_struct(F,F.haspng);
[tmp idx] = min(F.nseg);

if ~isempty(P.force_use_flowcell)
  fidx = grep(P.force_use_flowcell,F.segname,1);
  if length(fidx)>1, error('force_use_flowcell: multiple matches'); end
  if length(fidx)<1, error('force_use_flowcell: no match'); end
  fprintf('force_use_flowcell: %s\n',P.force_use_flowcell);
  idx = fidx;
end


fprintf('Linking to SegSeq output...\n');
force_link(F.segname{idx},[lawdir '/SegSeq_results.seg.txt']);
force_link(F.pngname{idx},[lawdir '/SegSeq_results.png']);

