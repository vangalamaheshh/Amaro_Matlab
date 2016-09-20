function [allpermst allpermsc] = exactclassperms(count,Nt,Nc,classes,permType)
%EXACTCLASSPERMS - generate random, class-constrained permutation indices

% process optional arguments
if ~exist('permType','var') || isempty(permType) || ~strcmp(permType,'trinomial')
    binomial = true;
end

Nclasses = length(classes);
permt = cell(1,Nclasses); % stores treatment permutation lists for each class
permc = cell(1,Nclasses); % stores control permutation lists for each class

% generate and store all permutations per-class in cell arrays
for c = 1:Nclasses
    if binomial
        [permt{c} permc{c}] = binomiter(Nt(c)+Nc(c),Nt(c));
    else
        [permt{c} permc{c}] = trinomiter(N(c),Nt(c),Nc(c));
    end
end
% The 'allperms' matrices contain the treatment and control parts of 
% the permutations, one permutation per column. 
allpermst = zeros(sum(Nt),count);
allpermsc = zeros(sum(Nc),count);
% iterate through the class permutaions odometer style,
% storing the concatenated permutations
codometer = ones(1,Nclasses);
digit = 1;
px = 1; % permutation index
while digit <= length(codometer)
    % concatenate current permutation column vector from class parts
    pt = []; pc = [];
    for c = 1:Nclasses
        if ~isempty(permc{c})
            pc = [pc;classes{c}( permc{c}(:,codometer(c)) )];
        end    
        if ~isempty(permt{c})
            pt = [pt;classes{c}( permt{c}(:,codometer(c)) )];
        end    
    end
    allpermst(:,px) = pt;
    allpermsc(:,px) = pc;
    px = px + 1;
    % increment class odometer with rollover
    digit = 1;
    while digit <= length(codometer) && codometer(digit) == size(permt{digit},2)
        codometer(digit) = 1;
        digit = digit + 1;
    end
    if digit <= length(codometer)
        codometer(digit) = codometer(digit) + 1;
    end
end

%%%%% SUBFUNCTIONS %%%%%
    
%% trinomial exact iterator function
function [perm1 perm2] = trinomiter(N,k1,k2)
% Return all ways of selecting 2 groups of k1 and k2 elements of
% 1:N. Returns two ki by N!/k1!k2!(N-k1-k2)! arrays (i=1,2) with the
% permutations in columns.
    
S = binomiter(N,k1+k2);     % all choices of combined group
T = binomiter(k1+k2,k1);    % all choices of subgroup

Nperm = size(S,2) * size(T,2);
perm1 = zeros(k1,Nperm);    % output space for the k1-sized group
perm2 = zeros(k2,Nperm);    % output space for the k2-sized group
for k = 1:size(T,2)
    % Trow1,Trow2 partition 1:k1+k2
    Trow1 = T(:,k);         % iterate over each subgroup combination
    Trow2 = (1:k1+k2)';
    Trow2(Trow1) = [];      % (fast difference set)
    if ~isempty(Trow1)
        perm1(:,(k-1)*size(S,2)+(1:size(S,2))) = S(Trow1,:);
    end
    if ~isempty(Trow2)
        perm2(:,(k-1)*size(S,2)+(1:size(S,2))) = S(Trow2,:);
    end
end

%% binomial exact iterator function
function [perms rest] = binomiter(N,k)
% Return all ways of selecting a group of k elements of 1:N
% result is a k by N!/k!(N-k)! array with the permutations
% in columns. Like nchoosek(1:N,k), but sideways and way faster.
if k == 0
    perms = zeros(0,1);
    rest = (1:N)';
elseif k > N
    perms = zeros(k,0);
    rest = [];
else
    % initialize output array(s) and its column iterator
    nck = nchoosek(N,k);
    perms = zeros(k,nck);
    if nargout > 1
        rest = zeros(N-k,nck);
    end
    ip = 1;

    % initialize iterator array
    ia = (1:k)';

    while true
        % output initial permutation
        perms(:,ip) = ia;
        if nargout > 1 && ~isempty(rest)
            rest(:,ip) = setdiff((1:N)',ia);
        end
        ip = ip + 1;
        % iterate to next permutation
        iap = k;
        ia(k) = ia(k)+1;
        while ia(iap) > N + iap - k
            iap = iap - 1;
            if iap == 0
                return
            else
                ia(iap:end) = ia(iap)+(1:(k-iap+1));
            end
        end
    end
end

