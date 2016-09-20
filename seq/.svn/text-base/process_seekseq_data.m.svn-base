function X = process_seekseq_data(instem)

T = load_seekseq_data([instem '-Tumor']);
N = load_seekseq_data([instem '-Normal']);

cand = T.W2R;
cand(:,5) = nan;
matchmarg = 6000;
for i=1:size(cand,1)
  idx = find(N.W2R(:,1)==cand(i,1) & N.W2R(:,2)<=cand(i,3)+matchmarg & N.W2R(:,3)>=cand(i,2)-matchmarg);
  if ~isempty(idx), cand(i,5) = idx(1); end
end
cand = cand(isnan(cand(:,5)),:);
[tmp ord] = sort(cand(:,4),'descend');
cand = cand(ord,:);

X.T=T;
X.N=N;
X.cand = cand;

Y=[];
Y.chr = cand(:,1);
Y.min = cand(:,2);
Y.max = cand(:,3);
Y.span = cand(:,3)-cand(:,2)+1;
Y.tumreads = cand(:,4);
save_struct(Y,[instem '.candidates.txt']);
