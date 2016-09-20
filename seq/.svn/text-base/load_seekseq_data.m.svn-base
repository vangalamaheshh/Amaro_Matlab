function X = load_seekseq_data(stem)

R = cell(24,1); W = cell(24,1);
for c=1:24, fprintf('chr%d ',c');
  R{c} = load_matrix([stem '.chr' num2str(c) '.run.txt']);
  W{c} = load_matrix([stem '.chr' num2str(c) '.weird.txt']);
end,fprintf('\n');
R = cat(1,R{:});
W = cat(1,W{:});
W = W(:,[1 3 2 4 6 5 7 8]);
W2 = sortrows(W(:,[4 5 6 6 7 8 1 2 3]));  % (with end2 first)
W2(:,4) = nan;  % map weird-pairs to the runs found in the normal
for i=1:size(W2,1)
  idx = find(R(:,1)==W2(i,1) & R(:,2)<=W2(i,2) & R(:,3)>=W2(i,2));
  if ~isempty(idx), W2(i,4) = idx; end
end

% find runs of weird pairs whose pairmates don't map to any run of normal pairs
thresh = 1000;
chr = -1;
first = -1;
last = -1;
n = 0;
nreg = 0;
W2R = [];

for i=1:size(W2,1)
  if chr==-1
    chr = W2(i,1);
    first = W2(i,2);
    n = 0;
    nreg = 0;
  else
    if W2(i,1)~=chr || W2(i,2)>last+thresh 
      W2R = [W2R;chr first last n nreg];
      chr = -1;
    end
  end    
  last = W2(i,2);
  if ~isnan(W2(i,4)), nreg=nreg+1; end
  n=n+1;
end

max_nreg = 0;
min_n = 2;
W2R = W2R(W2R(:,5)<=max_nreg & W2R(:,4)>=min_n,1:4);   % keep only the ones that dont map to any run of normal pairs

X=[];
X.R=R;
X.W=W;
X.W2=W2;
X.W2R=W2R;


