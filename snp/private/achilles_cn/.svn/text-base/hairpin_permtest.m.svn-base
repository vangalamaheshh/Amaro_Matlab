function R = hairpin_permtest(statfun,gscores,gtreat,gcontrol,classes,options)
%HAIRPIN_PERMTEST perform permutation tests on a set of unstructured hairpin data
%
% Perform permutation tests on a set of hairpin scores using differential
% statistics between two groups. The groups may, in general, be defined as
% different for each gene grouping of hairpins. The permutations may be
% constrained to only randomize scores within sample lineages.
%
%	R = HAIRPIN_PERMTEST(STATFUN,GSCORES,GTREAT,GCONTROL,CLASSES,OPTIONS)
%
% STATFUN is a statistical function used to evaluate the scores observed
% between the two groups defined by GTREAT and GCONTROL. CLASSES defines
% the lineage restriction, and OPTIONS holds optional arguments described
% with more detail below. The returned structure R is a structure with
% fields containing vectors of statistical results, one element per hairpin:
%
%   R.stat are the results of statfun on the observed scores.
%   R.ntreat and R.ncontrol are the number of samples in the treatment and
% control groups respectively.
%   R.mean_t and R.mean_c are the means of the scores in the treatment and
% control groups.
%   R.nperms are the number of permutations possible given the structure of
% the data defined by GTREAT, GCONTROL, and CLASSES.
%   R.pvalue and R.qvalue are the p- and q-values obtained from the
% permutation tests. In the "leftright" mode of operation, left- and right-
% tailed values are calculated from the same permutations, and returned
% instead in the fields R.pvalueL, R.pvalueR, R.qvalueL and R.qvalueR
%	R.type returns a character code for each p-value indicating how it was
% arrived at: '=' means all permutations were tested and the result is
% exact; '~' means Monte Carlo permutations were used and the result is
% approximate; '<' also means Monte Carlo, but no permutations surpassed
% the observed and a pseudocount was added; ' ' (space) means no
% permutation test was performed because the group was too small. When no
% permutations are performed for a hairpin, the p- and q-values are NaNs.
%
% The STATFUN inpu argument is the handle of a function that will generate
% statistics for the  permuted scores, e.g. difference of the mean. The
% prototypical STATFUN takes two 2- or two 3-dimensional array arguments
% for treatment and control classes. The first dimension enumerates scores
% within the treatment and control groups, and is generally different in
% size between the two input arrays. The second dimension indexes the
% permutations. The statistical function should compare the two arrays and
% reduce the first dimension using the statistic returning a row vector the
% same size as the second dimension in the 2-dimensional case. The
% 3-dimensional case is for parallelizing multiple hairpins with the same
% GTREAT/GCONTROL (gene) structure, the result should be a 3-dimensional
% array with a first dimension size of 1 and the 2nd/3rd dimension sizes
% matching the 2nd/3rd dimension sizes of the inputs. 
%
% GSCORES is an MxN data array of hairpin scores for M hairpins and N
% samples. GTREAT and GCONTROL are corresponding MxN arrays of logical
% values indicating whether each score belongs to the treatment or control
% group. A score cannot belong to both groups, however it can belong to
% neither ("no call"). 
%
% CLASSES is a cell array vector that partions the sample columns of the
% input into lineages (or other classes) that are permuted independently
% from one another. Each element contains a vector of indices of samples
% belonging to each lineage. Each sample column should be in exactly one
% of the elements of CLASSES.
%
% OPTIONS is a structure containing optional arguments.
%
%   OPTIONS.min_group defines the minimum number of samples required to be
% in a group before permutations are run. The default is 4.
%   OPTIONS.max_exact specifies the threshold number of permutations
% where the method switches from doing all permutations exhaustively to 
% doing a Monte Carlo sampling of permutations. The default is 10000.
%   OPTIONS.nTrials specifies how many random permutations are sampled if
% the Monte Carlo method is used. The default is 50000.
%   OPTIONS.tail specifies the directionality of the test: 'left' tests the
% fraction that control samples have a lower or equal statistic than the
% treatment sample; 'right' tests higher or equal; and 'both' (the default)
% returns the lower of 'left' or 'right' and but doubles the p-value;
% 'leftright' returns both values so that two hypotheses may be tested with
% the same permutations. OPTIONS.max_chunk is a ceiling for the number of
% hairpins with the same pattern GTREAT/GCONTROL group pattern used to limit
% memory usage by the function. The default is 100.
%   OPTIONS.permType allows the 'trinomial' variation of the permutation
% test where the "no-call" samples are permuted with the treatment and
% control samples. The default 'binomial' mode disregards no-calls.

    %% supply default options
    if ~exist('options','var') || isempty(options)
        options = struct;
    end
    if ~isfield(options,'min_group') || isempty(options.min_group)
        options.min_group = 4;
    end
    if ~isfield(options,'max_chunk') || isempty(options.max_chunk)
        options.max_chunk = 100;
    end

    %% calculate some statistics
    nhp = size(gscores,1);
    ntreat = sum(gtreat,2);
    ncontrol = sum(gcontrol,2);
    % calculate means (belong here???)
    mean_t = nan(nhp,1);
    mean_c = nan(nhp,1);
    for k = 1:size(gscores,1)
        mean_t(k) = nanmean(gscores(k,gtreat(k,:)));
        mean_c(k) = nanmean(gscores(k,gcontrol(k,:)));
    end
    
    %% initialize permutation input array
    [~,uid,gid] = unique([gtreat,gcontrol],'rows');
    premap = cell(1,length(uid));
    for k=1:length(uid)
        premap{k} = find(gid==k);
    end
    
    sizes = cellfun(@(x) size(x,1),premap);
    chunks = ceil(sizes/options.max_chunk);
    chunks(ntreat(uid)<options.min_group|ncontrol(uid)<options.min_group) = 0;
    
    nchunks = sum(chunks);
    slices  = cell(1,nchunks);
    putback = cell(1,nchunks);
    nanout = cell(1,nchunks);
    charout = cell(1,nchunks);
    treat = cell(1,nchunks);
    control = cell(1,nchunks);
    chunki = 1;
    for k = 1:length(chunks)
        if chunks(k) > 0
            slicex = premap{k};
            while ~isempty(slicex)
                % cut slicey from slicex
                ii = 1:min(length(slicex),options.max_chunk);
                slicey = slicex(ii);
                slicex(ii) = [];
                % set values of input and output fields
                treat{chunki} = gtreat(uid(k),:);
                control{chunki} = gcontrol(uid(k),:);
                slices{chunki} = gscores(slicey,:);
                nanout{chunki} = nan(length(slicey),1);
                charout{chunki} = repmat(' ',length(slicey),1);
                putback{chunki} = slicey;
                % next chunk
                chunki = chunki + 1;
            end
        end
    end
    seeds = num2cell(randi(1000,1,length(slices)));
    % build permutation input and output struct arrays 
    permin = struct('scores',slices,'treat',treat,'control',control,'seed',seeds);
    permout = struct('leftP',nanout,'rightP',nanout,'count',nanout,'stat',nanout,'type',charout);
    
        
    %% do the permutations
    %! (to permute in parallel, uncomment the matlabpool statements
    %! and change the 'for' to a 'parfor')
    %!matlabpool open 8
    % loop across matched genes
    Ntests = 0;
    for g = 1:length(permin)
        rand(permin(g).seed,1);
        verbose('Processing permutation %d of %d (%d hairpins).\n',30,g,length(permin),size(permin(g).scores,1));
        tmap = permin(g).treat;   % logical map of current treatment group
        cmap = permin(g).control; % logical map of current control group

        [permout(g).leftP,permout(g).count,permout(g).rightP,permout(g).type,permout(g).stat] = ...
            permtest(statfun,permin(g).scores,tmap,cmap,classes,options);
        Ntests = Ntests+1;
    end
    %!matlabpool close
    verbose('Performed %d permutation tests.',20,Ntests);

%% output results in input order
    R = struct;
  
    if strcmp(options.tail,'leftright')
        R.pvalueL = nan(nhp,1);
        R.pvalueR = nan(nhp,1);
    else
        R.pvalue = nan(nhp,1);
    end
    R.type = repmat(' ',nhp,1);
    R.stat = nan(nhp,1);
    % put values from slices back into linear array
    for k=1:length(putback)
        slicex = putback{k};
        if strcmp(options.tail,'leftright')
            R.pvalueL(slicex) = permout(k).leftP;
            R.pvalueR(slicex) = permout(k).rightP;
        else
            R.pvalue(slicex) = permout(k).leftP; % note: may be right if options.tail = 'right'
        end
        R.nperms(slicex) = permout(k).count;
        R.stat(slicex) = permout(k).stat;
        R.type(slicex) = permout(k).type;
    end
    R.ntreat = ntreat;
    R.ncontrol = ncontrol;
    R.mean_t = mean_t;
    R.mean_c = mean_c;
    
    % calculate false discovery rate
    if strcmp(options.tail,'leftright')
        permed = ~isnan(R.pvalueL);
        R.qvalueL = R.pvalueL;
        R.qvalueL(permed) = calc_fdr_value(R.pvalueL(permed));
        permed = ~isnan(R.pvalueR);
        R.qvalueR = R.pvalueR;
        R.qvalueR(permed) = calc_fdr_value(R.pvalueR(permed));
    else
        permed = ~isnan(R.pvalue);
        R.qvalue = R.pvalue;
        R.qvalue(permed) = calc_fdr_value(R.pvalue(permed));
    end
    
