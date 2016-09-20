function fle2 = find_pairmate(fle1)

lanes = load_lanelist;

[f1 fi fj] = unique(fle1);

f2 = nan(size(f1));
for i=1:length(f1)
  if f1(i)<1 || f1(i)>slength(lanes), fprintf('%d out of range\n',f1(i));
  else
    idx = find(strcmp(lanes.FC,lanes.FC{f1(i)}) & ...
      strcmp(lanes.lane,lanes.lane{f1(i)}) & ...
      ~strcmp(lanes.end,lanes.end{f1(i)}));
    if ~isempty(idx), f2(i) = lanes.FLE(idx);
    else fprintf('%d not found\n',f1(i));
end,end,end

fle2 = f2(fj);
