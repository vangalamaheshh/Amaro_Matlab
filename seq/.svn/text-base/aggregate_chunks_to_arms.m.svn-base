function Xarm = aggregate_chunks_to_arms(X,P)
% Mike Lawrence 2009-07-18

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'omit_arms',{});

P.omit_arms = regexprep(P.omit_arms','X','23');
P.omit_arms = regexprep(P.omit_arms','Y','24');

cen = load_cen;
cen = mean(cen,2);

Xarm = [];
Xarm.lane = X.lane;

arm = 1;

for c = 1:24
  cidx = find(X.chunk.chr == c);
  pidx = cidx(find(X.chunk.end(cidx)<cen(c)));
  qidx = setdiff(cidx,pidx);
  for a = 1:2
    if a==1, armstring = [num2str(c) 'p']; else armstring = [num2str(c) 'q']; end
    if ismember(armstring,P.omit_arms), continue; end
    if a==1, idx=pidx; else idx=qidx; end
    Xarm.chunk.num(arm,1) = arm;
    Xarm.chunk.chr(arm,1) = c;
    Xarm.chunk.start(arm,1) = min(X.chunk.start(idx));
    Xarm.chunk.end(arm,1) = max(X.chunk.end(idx));
    Xarm.dat(arm,:) = sum(X.dat(idx,:));
    arm=arm+1;
  end
end

[Xarm.nchunks, Xarm.nlanes] = size(Xarm.dat);

