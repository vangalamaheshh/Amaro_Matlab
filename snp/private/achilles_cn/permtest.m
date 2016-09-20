function [pval count oval type stat] = permtest(statfun,scores,gtreat,gcontrol,classes,options)
%PERMTEST Perform a statistical permutation test on data.
%
% [PVAL COUNT OVAL TYPE STAT] = permtest(STATFUN,SCORES,GTREAT,GCONTROL,CLASSES,OPTIONS);
%
% Returns p-values in PVALs using the statistical function handle STATFUN 
% on two groups of data within the row vector SCORES calculated by
% permuting how the samples are partitioned into the two groups. The
% groups are designated by the logical row vectors GTREAT and
% GCONTROL which are true where the corresponding data element is
% in the group. Depending on the OPTIONS, the permutation test may be a
% Monte Carlo test or an exhaustive exact test of all permutations.
%
% The optional COUNT output returns the theoretical number of permutations
% possible given the structure of the data defined by GTREAT, GCONTROL, and
% CLASSES. OVAL returns the opposite p-value from PVAL, that is, if OPTIONS
% specify right p-values, OVAL are the left p-values and vice-versa. It is
% a copy of PVAL in the two-tailed case. TYPE returns a character for each
% p-value indicating how it was arrived at: '=' means all permutations
% were used and the result is exact; '~' means Monte Carlo permutations
% were used and the result is approximate; '<' also means Monte Carlo, but
% no permutations surpassed the observed and a pseudocount was added. STAT
% is the observed statistic for each p-value.
%
% STATFUN is the handle of a function that will generate statistics for the 
% permuted scores, e.g. difference of the mean. The prototypical STATFUN
% takes two 2- or two 3-dimensional array arguments for treatment and
% control classes. The first dimension enumerates scores within the
% treatment and control groups, and is generally different in size between
% the two input arrays. The second dimension indexes the permutations. The
% statistical function should compare the two arrays and reduce the first
% dimension using the statistic returning a row vector the same size as the
% second dimension in the 2-dimensional case. The 3-dimensional case is for
% parallelizing multiple hairpins with the same GTREAT/GCONTROL (gene)
% structure, the result should be a 3-dimensional array with a first
% dimension size of 1 and the 2nd/3rd dimension sizes matching the 2nd/3rd
% dimension sizes of the inputs. 
%
% The optional CLASSES input allows the samples to potentially be
% partitioned into several lineage classes within which permutations are
% performed independently. Each cell in CLASSES is a set of indices for
% samples of the same lineage. If CLASSES is empty, lineage restrictions
% are not observed.
%
% The optional OPTIONS structure specifies details of the permutation test.
% The OPTIONS.max_exact field specifies the threshold number of
% permutations where the method switches from doing all permutations to
% doing a Monte Carlo sampling of permutations. The OPTIONS.nTrials field
% specifies how many random permutations are sampled if the Monte Carlo
% method is used. The OPTIONS.tail specifies the directionality of the
% test: 'left' tests the fraction that control samples have a lower or
% equal statistic than the treatment sample; 'right' tests higher
% or equal; and 'both' returns the lower of 'left' or 'right', but
% doubles the p-value.
    
    % optional argument handling
    if ~exist('classes','var') || isempty(classes)
        classes = {1:size(scores,2)};
    end
    if ~exist('options','var') || isempty(options)
        options = struct;
    end
    if ~isfield(options,'tail') || isempty(options.tail)
        options.tail = 'both';
    end
    if ~isfield(options,'max_exact') || isempty(options.max_exact)
        options.max_exact = 10000;
    end
    if ~isfield(options,'nTrials') || isempty(options.max_exact)
        options.nTrials = 50000;
    end
    if ~isfield(options,'permType') || isempty(options.permType)
        options.permType = 'binomial';
    end

    binomial = strcmp(options.permType,'binomial'); % TODO catch bad values
    
    classes = cellfun(@(x) x(:),classes,'UniformOutput',false);
    
    %% for binomial statistics, reduce input to just treatment and control elements
    if binomial
        reducer = gtreat|gcontrol;
        scores = scores(:,reducer);
        gtreat = gtreat(reducer);
        gcontrol = gcontrol(reducer);
        % reduce classes partition description
        redmap = cumsum(reducer)';
        classes = cellfun(@(c) subsref(redmap(c),substruct('()',{reducer(c)})),...
                    classes,'UniformOutput',false);
    end
    
    %% calculate statistic for the actual obtained result
    % first, put indices in same order that permutations will be in to
    % mitigate roundoff differences
    permorder = cell2mat(classes');
    pscores = scores(:,permorder);
    real_statval = feval(statfun,pscores(:,gtreat(permorder))',pscores(:,gcontrol(permorder))');
    
    %% count permutations for each class
    [count N Nt Nc] = countclassperms(classes,gtreat,gcontrol,binomial);
    
    %% generate permutation indices with class constraints
    nhp = size(scores,1);
    if count > options.max_exact
        % Monte Carlo permutations
        [allpermst allpermsc] = randclassperms(N,Nt,Nc,classes,options.nTrials,options.permType);
        type = repmat('~',nhp,1);
    else
        % exact permutations 
        [allpermst allpermsc] = exactclassperms(count,Nt,Nc,classes,options.permType);
        type = repmat('=',nhp,1);
    end
    
    %% calculate statistics for all permutations and derive p-value
    % permute scores in treatment and control groups
    t_scores = scores(:,allpermst);
    c_scores = scores(:,allpermsc);
    % rearrange so that the optional first dimension of the scores 
    % forms the optional last dimension of resulting 3-dimensional array.
    t_scores = permute(reshape(t_scores,[nhp,size(allpermst)]),[2,3,1]);
    c_scores = permute(reshape(c_scores,[nhp,size(allpermsc)]),[2,3,1]);

%! as of 10/21, any filtering will be done by the passed in statistical
%! function by returning a statistic of NaN
%!    t_means = mean(t_scores,1);
%!    keepers = t_means >= options.lowerMeanLimit ...
%!            & t_means <= options.upperMeanLimit;

    % apply statistical function
    perm_statvals = reshape(feval(statfun,t_scores,c_scores),size(allpermst,2),nhp);
    % sum
    pLeft = nan(nhp,1);
    pRight = nan(nhp,1);
    for h = 1:nhp
        pLeft(h) = nansum(perm_statvals(:,h) <= real_statval(h),1);
        pRight(h) = nansum(perm_statvals(:,h) >= real_statval(h),1);
    end
    % pseudocounts
    type(pLeft==0|pRight==0) = '<';
    pLeft(pLeft==0) = 1;
    pRight(pRight==0) = 1;
    % frequency
    nGoodPerms = sum(~isnan(perm_statvals))';
    pLeft = pLeft ./ nGoodPerms;
    pRight = pRight ./ nGoodPerms;
    
    % return pvalues, other values
    if strcmp(options.tail,'left') || strcmp(options.tail,'leftright')
        pval = pLeft;
        oval = pRight;
    elseif strcmp(options.tail,'right') || strcmp(options.tail,'rightleft')
        pval = pRight;
        oval = pLeft;
    else % default 'both'
        pval = 2 * pLeft;
        pval(pRight < pLeft) = pRight;
        oval = pval;
    end
    stat = real_statval;

