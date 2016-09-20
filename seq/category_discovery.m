% CONVENIENT STARTING FORMAT:  16x6 matrix

% ROWS:      COLS:  A_A A_C A_G A_T | C_A ... | ... | ... T_G T_T
% A->C
% A->G
% A->T
% C->A
% C->G
% C->T

% TSP data ("pooled" from silent + singleton nonsilent, and correctly strand-collapsed)

TSP.n = [0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1;0,1,1,2,0,3,3,2,1,1,1,1,0,1,1,4;1,2,0,1,2,1,5,3,0,1,0,0,0,1,1,0;1,4,4,1,5,6,10, ...
     6,5,7,6,2,2,9,5,3;0,1,0,0,2,1,1,2,0,1,3,1,4,4,2,4;10,4,4,4,9,11,8,10,3,6,6,1,13,9,14,7];
TSP.N = [11005,7627,13022,7888,11656,12699,22157,11737,12887,9198,14937,9276,5574,5916,5861,5710;12851,8868,14665,9328, ...
     12578,13409,21675,13574,13551,9618,14757,10113,5978,6221,5697,6464;11170,7627,12697,7888,12062,12700,21675,11737,13132, ...
          9198,14503,9276,5762,5916,5702,5710;11115,9436,4769,8651,19393,14734,9097,14893,14413,14766,7208,13627,14294,14685, ...
          5454,13093;11115,9436,4612,8651,19393,14734,8631,14893,14413,14765,6830,13626,13819,14297,5043,12843;12935, ...
          10999,5380,11091,19702,14935,8631,18140,15306,15636,7127,16869,14926,15288,5383,15694];
TSP.N = TSP.N * 188;

% TCGA data ("pooled" from silent + non-top-40 nonsilent, and correctly strand-collapsed)

TCGA.n = [0,0,0,0,2,1,0,0,0,0,0,0,0,0,1,0;3,0,1,1,0,0,2,2,1,1,0,0,0,1,1,0;0,0,2,0,0,0,0,0,1,0,0,0,0,0,0,0;0,1,1,1,2,1,0,1, ...
                  1,1,4,0,2,1,0,0;0,1,0,0,1,2,0,4,1,1,0,2,0,1,1,2;0,6,16,2,7,6,15,8,2,8,16,6,2,4,5,11];
TCGA.N = [6755,3538,15125,2908,12556,12763,32975,8697,10820,7491,20145,4938,4878,5881,7119,4841;20151,13447,20482, ...
          14995,14996,16703,27379,16052,16737,11819,17991,12773,10102,7748,5778,7578;8027,3538,9333,2908,13567,12763,27379, ...
          8697,11849,7491,15350,4938,5662,5881,5795,4841;6077,3811,4980,2666,16184,12856,10922,11418,9406,10605,7686, ...
          5958,10618, 10798,7024,7943;6077,3811,3061,2666,16184,12856,8733,11418,9406,10605,5822,5958,9183,9638,5202,6762; ...
          13085,11134,5378,14774,17424,14300,8733,21904,14498,16614,7212,21072,22424,14824,6027,20305];
TCGA.N = TCGA.N * 87;

% 2009-12-05 MELANOMA 20xC2K from: /xchip/tcga_scratch/lawrence/mel/analysis/20091205/run.m
Y = load_struct('/xchip/tcga_scratch/lawrence/mel/analysis/20091205/mut_categs.txt','%f%s%f%f%f%f%f');
MEL.n = [Y.C(1:16)';Y.G(1:16)';Y.T(1:16)';Y.A(17:32)';Y.G(17:32)';Y.T(17:32)'];
MEL.N = [repmat(Y.N(1:16)',3,1);repmat(Y.N(17:32)',3,1)]

% CHOOSE WHICH DATA SET TO USE

%orig_n = TSP.n; orig_N = TSP.N;
%orig_n = TCGA.n; orig_N = TCGA.N;
orig_n = MEL.n; orig_N = MEL.N;

% reformat into 4-dimensional matrix

% dimensions:
%   (1)  5' base (1234=ACGT)
%   (2)  3' base (1234=ACGT)
%   (3)  "from" base (12=AC)
%   (4)  "to" outcome (1=transition, 2=flip_transversion, 3=skew_transversion)

n = zeros(4,4,2,3);
N = zeros(4,4,2,3);

for base5 = 1:4
  for base3 = 1:4
    for oldbase = 1:2
      for muttype = 1:3
        col = 4*(base5-1)+base3;
        switch oldbase
          case 1 % A
            switch muttype
              case 1, row = 2;  % A->G transition
              case 2, row = 3;  % A->T transversion(flip)
              case 3, row = 1;  % A->C transversion(skew)
            end
          case 2 % C
            switch muttype
              case 1, row = 6;  % C->T transition
              case 2, row = 5;  % C->G transversion(flip)
              case 3, row = 4;  % C->A transversion(skew)
            end
        end
        n(base5,base3,oldbase,muttype) = orig_n(row,col);
        N(base5,base3,oldbase,muttype) = orig_N(row,col);
      end
    end
  end
end

% test split

unsplit = {[1:4],[1:4],[1:2],[1:3]};
at = {[1:4],[1:4],[1],[1:3]};
cg = {[1:4],[1:4],[2],[1:3]};

initial = {unsplit};
split1 = {at;cg};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% navigate breakdown tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial = {unsplit};
fprintf('INITIAL\n\n');
H_initial = entropy_by_parts(n,N,initial);
stats = reportrule(n,N,initial);
std_tot = sum(stats.ci(:,2)-stats.rate)/1.98;
fprintf('\n    H_initial = %d  std_tot = %d\n\n', H_initial, std_tot);

current = initial;
step=1;
while(1)
  fprintf('STEP %d\n\n',step);
  H_current = entropy_by_parts(n,N,current);
  best_new = {};
  best_dH = 0;
  for t=1:length(current)
    rest = current(setdiff(1:length(current),t));
    parent = current{t};
    for d=1:4
      tobreak = parent{d};
      p = powerset(tobreak);
      minel = min(tobreak);
      for i=2:length(p)-1
        if ismember(minel,p{i})
          child1 = parent;
          child2 = parent;
          child1{d} = p{i};
          child2{d} = setdiff(tobreak, p{i});
          new_parts = [rest;{child1};{child2}];
          H_new = entropy_by_parts(n,N,new_parts);
          dH = H_new - H_current;
          if dH<best_dH
            best_dH = dH;
            best_new = new_parts;
  end,end,end,end,end

  stats = reportrule(n,N,best_new);
  std_tot = sum(stats.ci(:,2)-stats.rate)/1.98;
  fprintf('\n    dH = %d  H_new = %d  std_tot = %d\n\n', best_dH, H_current+best_dH, std_tot);
  current = best_new;
  step=step+1;
  if step>20, break; end
end % next step


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find best possible category set of a given size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial = {unsplit};
leaves = {initial};

for size=1:10
  fprintf('k=%d:  ', size);

  if size>1
    old_leaves = leaves;
    leaves = {};
    fprintf(' [gen] ');
    for l = 1:length(old_leaves)
      leaf = old_leaves{l};
      for t=1:length(leaf)
        rest = leaf(setdiff(1:length(leaf),t));
        parent = leaf{t};
        for d=1:4
          tobreak = parent{d};
          p = powerset(tobreak);
          minel = min(tobreak);
          for i=2:length(p)-1
            if ismember(minel,p{i})
              child1 = parent;
              child2 = parent;
              child1{d} = p{i};
              child2{d} = setdiff(tobreak, p{i});
              new_leaf = [rest;{child1};{child2}];
              leaves = [leaves;{new_leaf}];
    end,end,end,end,end
    leaves = remove_duplicate_leaves(leaves);
  end

  fprintf('%d possible category set', length(leaves));
  if length(leaves)>1,fprintf('s');end
  fprintf('\n');

  fprintf(' [pickbest] ');
  best_l = 0;
  best_H = Inf;
  for l=1:length(leaves)
    H = entropy_by_parts(n,N,leaves{l});
    if H<best_H
      best_l = l;
      best_H = H;
    end
  end

  % report results
  fprintf('Best:\n');
  best_leaf = leaves{best_l};
  stats = reportrule(n,N,best_leaf);
  std_tot = sum(stats.ci(:,2)-stats.rate)/1.98;
  fprintf('\n    H_final = %d  std_tot = %d\n\n', best_H, std_tot);
end % next size


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find best possible category set of a given size
% OPTIMIZED VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial = {unsplit};
leaves = {initial};

estimated_growth_factor = 30;   % actual split closer to 20, but let's err on generous side
for size=1:10
  fprintf('k=%d:  ', size);
  if size>1
    old_leaves = leaves;
    typical_leaf = repmat({unsplit},size,1);
    leaves = repmat(typical_leaf,length(old_leaves)*estimated_growth_factor,1);
    lastleaf = 0;
    fprintf(' [gen] ');
    for l = 1:length(old_leaves)
      if ~mod(l,10), fprintf(' [leaf %d/%d] ',l,length(old_leaves)); end
      leaf = old_leaves{l};
      for t=1:length(leaf)
        rest = leaf(setdiff(1:length(leaf),t));
        parent = leaf{t};
        for d=1:4
          tobreak = parent{d};
          p = powerset(tobreak);
          p = p(2:end-1);   % remove empty and full sets
          np = length(p);
          new_leaf = [rest;{parent};{parent}];
          new_leaves = repmat({new_leaf},np,1);
          child1i = length(new_leaf);
          child2i = child1i-1;
          for i=1:np
            new_leaves{i}{child1i}{d} = p{i};
            new_leaves{i}{child2i}{d} = setdiff(tobreak,p{i});
          end
          leaves(lastleaf+1:lastleaf+length(new_leaves)) = new_leaves;
          lastleaf = lastleaf + length(new_leaves);
    end,end,end
    leaves = remove_duplicate_leaves(leaves(1:lastleaf));
  end
  fprintf('%d possible category set', length(leaves));
  if length(leaves)>1,fprintf('s');end
  fprintf('\n');

  fprintf(' [pickbest] ');
  best_l = 0; best_H = Inf;
  for l=1:length(leaves)
    if ~mod(l,1000), fprintf(' %d/%d ',l,length(leaves)); end
    H = entropy_by_parts(n,N,leaves{l});
    if H<best_H, best_l = l; best_H = H; end
  end

  % report results
  fprintf('Best:\n');
  best_leaf = leaves{best_l};
  stats = reportrule(n,N,best_leaf);
  std_tot = sum(stats.ci(:,2)-stats.rate)/1.98;
  fprintf('\n    H_final = %d  std_tot = %d\n\n', best_H, std_tot);
end % next size



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new version based on all-integer processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unsplit = categ_to_int({[1 2 3 4] [1 2 3 4] [1 2] [1 2 3]});
mask = uint16(15*16.^[0:3]);
nybs = uint16((1:14)'*(16.^[0:3]));

initial = unsplit;
leaves = initial;
for k=1:10
  fprintf('k=%d\n',k); tic
  if k>1
    old_leaves = leaves;
    leaves = zeros(34^(k-1),k,'uint16');
    lastleaf = 0;
    for l=1:length(old_leaves)
      if ~mod(l,1000), fprintf(' [leaf %d/%d] ',l,length(old_leaves)); end
      leaf = old_leaves(l,:);
      for c=1:k-1   % choose which category to split
        parent = leaf(c);
        for d=1:4   % choose which dimension to split along
          tobreak = bitand(parent,mask(d));
          tokeep = parent-tobreak;
          frags1 = bitand(tobreak,nybs(:,d));
          frags2 = tobreak-frags1;
          frags = [frags1 frags2];
          frags(~frags1|~frags2,:)=[];
          children = tokeep+frags;
          nc = size(children,1);
          leaves(lastleaf+1:lastleaf+nc,1:k-1) = repmat(leaf,nc,1);
          leaves(lastleaf+1:lastleaf+nc,[c k]) = children;
          lastleaf = lastleaf + nc;
    end,end,end
    fprintf('before condensing: %d\n',lastleaf);
    leaves = unique(sort(leaves(1:lastleaf,:),2),'rows');
  end
  fprintf('%d possible category set', size(leaves,1));
  if length(leaves)>1,fprintf('s');end
  fprintf('\n'); toc
end

