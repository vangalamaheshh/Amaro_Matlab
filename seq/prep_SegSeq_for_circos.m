function prep_SegSeq_for_circos(infile,outfile,P)
if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'omit_chrX_cn',false);
P = impose_default_value(P,'omit_chrY_cn',false);
P = impose_default_value(P,'species','hs');

mincap = 0;
maxcap = 4;

% load Derek's seg file
X = load_SegSeq_results(infile);
if P.omit_chrY_cn
  X = reorder_struct(X,X.chr~=24);
end
if P.omit_chrX_cn
  X = reorder_struct(X,X.chr~=23);
end

% create circos heatmap data file
Y = cell(slength(X),1);
if ismember(P.species,'mm')
    X.chromosome=regexprep(X.chromosome,'chr','')
end
for i=1:slength(X)
  Y{i} = sprintf('%s%s %d %d %f',P.species,...
    X.chromosome{i}, X.start(i), X.end(i),...
    max(mincap,min(maxcap,X.copyratio(i))));
end

save_lines(Y,outfile);


