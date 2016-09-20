function R = perform_permutations(R,statfun,gscores,gtreat,gcontrol,classes,options)
%
% TODO: head documentation
%
%

    % supply default options
    if ~exist('options','var') || isempty(options)
        options = struct;
    end
    if ~isfield(options,'min_group') || isempty(options.min_group)
        options.min_group = 4;
    end

    % initialize result arrays
    R.nperms = nan(length(R.acgene),1);
    if strcmp(options.tail,'leftright')
        R.pvalueR = nan(length(R.acgene),1);
        R.pvalueL = nan(length(R.acgene),1);
    else
        R.pvalue = nan(length(R.acgene),1);
    end
    R.ntreat = zeros(length(R.acgene),1);
    R.ncontrol = zeros(length(R.acgene),1);

    % loop across matched genes
    Ntests = 0;
    for g = 1:length(R.acgene)

        if ~mod(g,100)
            verbose('...processing gene %d of %d...%d permutation test so far',...
                    30,g,length(R.acgene),Ntests); 
        end
        ag = R.acgene(g);      % current achilles gene index
        cg = R.clgene(g);      % current CLE CN gene index
        tmap = gtreat(cg,:);   % logical map of current treatment group
        cmap = gcontrol(cg,:); % logical map of current control group

        Ntreatment = sum(tmap);
        Ncontrol = sum(cmap);
        R.ntreat(g) = Ntreatment;
        R.ncontrol(g) = Ncontrol;
        R.mean_t(g) = mean(gscores(ag,tmap));
        R.mean_c(g) = mean(gscores(ag,cmap));
        if all( [Ntreatment Ncontrol] >= options.min_group )
            if strcmp(options.tail,'leftright')
                [R.pvalueL(g) R.nperms(g) R.pvalueR(g)] = ...
                    permtest(statfun,gscores(ag,:),tmap,cmap,classes,options);
                Ntests = Ntests+1;
            else
                [R.pvalue(g) R.nperms(g)] = ...
                    permtest(statfun,gscores(ag,:),tmap,cmap,classes,options);
                Ntests = Ntests+1;
            end
        end
    end
    verbose('Performed %d permutation tests.',20,Ntests);
    
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
    
