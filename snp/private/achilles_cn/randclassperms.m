function [allpermst allpermsc] = randclassperms(N,Nt,Nc,classes,nTrials,permType)
%RANDCLASSPERMS - generate random, class-constrained permutation indices

% process optional arguments
if ~exist('nTrials','var') || isempty(nTrials)
    nTrials = 50000;
end
if ~exist('permType','var') || isempty(permType)
    binomial = true;
else
    binomial = ~strcmp(permType,'trinomial');
end

% The 'allperms' matrices contain the treatment and control parts of 
% the permutations, one permutation per column. 
allpermst = zeros(sum(Nt),nTrials);
allpermsc = zeros(sum(Nc),nTrials);
kt = 1;
kc = 1;
for c = 1:length(classes)
    % generate and store random permutations for a class
    if binomial
        [permt permc] = binorander(Nt(c)+Nc(c),Nt(c),nTrials);
    else
        [permt permc] = trinorander(N(c),Nt(c),Nc(c),nTrials);
    end
    if Nt
        allpermst(kt:kt+Nt(c)-1,:) = classes{c}(permt);
        kt = kt + Nt(c);
    end
    if Nc
        allpermsc(kc:kc+Nc(c)-1,:) = classes{c}(permc);
        kc = kc + Nc(c);
    end
end        

%%%%% SUBFUNCTIONS %%%%%
    
%% subfunction: trinomial Monte Carlo chooser function
function [perm1 perm2] = trinorander(N,k1,k2,trials)
% Return all ways of selecting 2 groups of k1 and k2 elements of
% 1:N. Returns two ki by N!/k1!k2!(N-k1-k2)! arrays (i=1,2) with the
% permutations in columns.
    deck = (1:N)';
    perm1 = zeros(k1,trials);
    perm2 = zeros(k2,trials);
    for tx = 1:trials
        for i = 1:k1
            j = i+floor((N-i+1)*rand);
            t = deck(j);
            deck(j) = deck(i);
            deck(i) = t;
            perm1(i,tx) = t;
        end
        for i = 1:k2
            j = k1+i+floor((N-k1-i+1)*rand);
            t = deck(j);
            deck(j) = deck(k1+i);
            deck(k1+i) = t;
            perm2(i,tx) = t;
        end
    end
    
%% subfunction: binomial Monte Carlo chooser
function [perms rest] = binorander(N,k,trials)
% Return ntrials random selections of a group of k elements of 1:N.
% Result is a k by trials array with each permutation in a column.
    if k == 0
        perms = zeros(0,trials);
        rest = repmat((1:N)',1,trials);
    elseif k > N
        perms = zeros(k,0);
        rest = [];
    else
        perms = zeros(k,trials);
        if nargout > 1
            rest = zeros(N-k,trials);
        end
        deck = (1:N)';
        % loop over trials
        for tx = 1:trials
            % Knuth-shuffle the deck
            for i = N:-1:2
                j = ceil(i*rand);
                t = deck(i);
                deck(i) = deck(j);
                deck(j) = t;
            end
            % pick the first k
            perms(:,tx) = deck(1:k);
            if nargout > 1 && N > k
                rest(:,tx) = deck(k+1:N);
            end
        end
    end
