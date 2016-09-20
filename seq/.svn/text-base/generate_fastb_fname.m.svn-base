function fname=generate_fastb_fname(i,freeze)
if ~exist('freeze','var') || isempty(freeze), freeze=default_freeze; end

lanes = load_lanelist;
freezes = load_freezelist;

if freeze<0 || freeze>slength(freezes), error('No such freeze'); end

idx = find(lanes.FLE==i);
if length(idx)~=1
  fprintf('Lane index not found in lanelist!\n');
  fname = 'none';
else
  if strcmp(lanes.date{idx},'None')
    fprintf('This lane was single-endread-only, and it lacks that endread.\n');
    fname = 'none';
  else
    if freeze==0
      fname=[ '/seq/solexa/pipelineOutput/' lanes.date{idx} '_' lanes.FC{idx} '/' ...
        lanes.FC{idx} '.' lanes.lane{idx} '.fastb' ];
    else 
    if strcmp(lanes.tumornormal{idx},'tumor'), fname = freezes.tumorbasedir{freeze};
    else fname = freezes.normalbasedir{freeze}; end
    fname = [fname '/alignments' ...
      '/' lanes.FC{idx} '.' lanes.lane{idx} '/' lanes.date{idx} '/' ...
      lanes.FC{idx} '.' lanes.lane{idx} '.fastb' ];

%%      fname = ['/home/radon01/kiran/tmp/tcga/tcga-' lanes.tumornormal{idx} '-freeze' ...
%%        num2str(freeze) '-' freezes.date{freeze} '/alignments' ...
%%        '/' lanes.FC{idx} '.' lanes.lane{idx} '/' lanes.date{idx} '/' ...
%%         lanes.FC{idx} '.' lanes.lane{idx} '.fastb' ];
    end
  end
end
