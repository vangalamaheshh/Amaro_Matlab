function R = perform_perms_parallel(R,statfun,gscores,gtreat,gcontrol,classes,options)
%
% TODO: head documentation
%
%

    % supply default options
    if ~exist('options','var') || isempty(options)
        options = struct;
    end
    if ~isfield(options,'nprocesses') || isempty(options.nprocesses)
        options.nprocesses = 4;
    end
    if ~isfield(options,'min_group') || isempty(options.min_group)
        options.min_group = 4;
    end

    % initialize result arrays
    nperms = nan(length(R.acgene),1);
    if strcmp(options.tail,'leftright')
        lPvals = nan(length(R.acgene),1);
        rPvals = nan(length(R.acgene),1);
    else
        pvals = nan(length(R.acgene),1);
    end
    ntreat = zeros(length(R.acgene),1);
    ncontrol = zeros(length(R.acgene),1);
    
    mean_t = nan(length(R.acgene),1);
    mean_c = nan(length(R.acgene),1);
    
    %! premap scores and groupings prior to permutations
    ascores = gscores(R.acgene,:);
    atreat = gtreat(R.clgene,:);
    acontrol = gcontrol(R.clgene,:);
    
    min_group = options.min_group;
    tail = options.tail;
    
    matlabpool open 6
    % loop across matched genes
    Ntests = 0;
    parfor g = 1:length(R.acgene)

%        if ~mod(g,100)
%            verbose('...processing gene %d of %d...%d permutation tests so far',...
%                    30,g,length(R.acgene),Ntests); 
%        end
        tmap = atreat(g,:);   % logical map of current treatment group
        cmap = acontrol(g,:); % logical map of current control group

        Ntreatment = sum(tmap);
        Ncontrol = sum(cmap);
        ntreat(g) = Ntreatment;
        ncontrol(g) = Ncontrol;
        mean_t(g) = mean(ascores(g,tmap));
        mean_c(g) = mean(ascores(g,cmap));
        if all( [Ntreatment Ncontrol] >= min_group )
            if strcmp(tail,'leftright')
                [lPvals(g) nperms(g) rPvals(g)] = ...
                    permtest(statfun,ascores(g,:),tmap,cmap,classes,options);
            else
                [pvals(g) nperms(g)] = ...
                    permtest(statfun,ascores(g,:),tmap,cmap,classes,options);
            end
            Ntests = Ntests+1;
        end
    end
    verbose('Performed %d permutation tests.',20,Ntests);
    matlabpool close
    
    R.ntreat = ntreat;
    R.ncontrol = ncontrol;
    R.mean_t = mean_t;
    R.mean_c = mean_c;
    R.nperms = nperms;
    if strcmp(options.tail,'leftright')
        R.pvalueR = rPvals;
        R.pvalueL = lPvals;
    else
        R.pvalue = pvals;
    end
    
    % calculate False Discovery Rate
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
    
