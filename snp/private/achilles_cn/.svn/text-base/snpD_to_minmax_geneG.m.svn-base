function G = snpD_to_minmax_geneG(D,rg,rg_field,options)
% SNPD_TOGENEG Convert SNP copy number data to gene level data
%
%   G = snpD_to_geneG(D,RG,RG_FIELD,options)
%
%   Converts SNP-level copy number data in D to gene level copy data in G
%   using the reference genome in structure array RG. RG_FIELD is a string
%   specifying the field in RG to use as a gene identifier, default 'symb'.
%   OPTIONS is a struct whose fields are optional parameters:
%       OPTIONS.collapse_method can be 'all' (default), 'mean', 'median',
%               'min', 'max', or 'extreme'
%       OPTIONS.find_snps_type is the type argument passed to find_snps
%               0 for SNP markers strictly contained in the gene region
%               1 to allow flanking SNPS to give gene CN data
%               2 to use flanking SNPS only if none are strictly contained
%

    if ~exist('options','var')
        options = struct;
    end
%!    if ~isfield(options,'collapse_method')
%!        options.collapse_method = 'all';
%!    end
    if ~isfield(options,'find_snps_type')
        options.find_snps_type = 1;
    end
%!    if ~isfield(options,'nan_thresh')
%!        options.nan_thresh = 0;
%!    end
    if ~exist('rg_field','var') || isempty(rg_field)
        rg_field = 'symb';
    end

    % compress genes with same identifier in reference genome
    ids = {rg.(rg_field)};
    uqids = unique(ids);
    [gmatch,m1,m2]=match_string_sets_hash(ids,uqids);
    % m1 indexes ids
    % m2 corresponds ids, indexes uqids
    
    % create mapping from compresssed reference genome to snps
    usize = [length(uqids) 1];
    snps = cell(usize);    % cell array of snp positions
    rgindex = cell(usize); % cell array of reference genome indices
    chrn = zeros(usize);   % chromosome (NaN if ambiguous)
    g_start = Inf(usize);  % minimum gene start
    g_end = zeros(usize);  % maximum gene end
    snp_start = Inf(usize);% minimum flanking snp start
    snp_end = zeros(usize);% maximum flanking gene end
    
    % loop over unique gene IDs
    for ag = 1:length(uqids)
        if mod(ag,50)==0
            disp(ag);
        end
        matches = find(gmatch(:,ag));
        rgindex{ag} = matches;
        if ~isempty(matches)
            % loop over all matches to unique gene ID
            for k = matches'
                g_start(ag) = min(g_start(ag),rg(k).start);
                g_end(ag) = max(g_end(ag),rg(k).end);
                if chrn(ag) == 0
                    chrn(ag) = rg(k).chrn;
                elseif chrn(ag) ~= rg(k).chrn;
                    chrn(ag) = NaN;
                    verbose('Conflicting chromosomes for gene %s !',30, ...
                            rg(k).(rg_field));
                end
                if g_end(ag) > g_start(ag) && ~isnan(chrn(ag))
                    snps{ag} = union(snps{ag},...
                                     find_snps(D,chrn(ag),g_start(ag),g_end(ag),1));
                else
                    verbose('Bad gene data for %s !',30, ...
                            rg(k).(rg_field));
                end
            end
            % get positions of (possibly flanking) snps
            if ~isempty(snps{ag})
                snp_start(ag) = D.pos(min(snps{ag}));
                snp_end(ag) = D.pos(max(snps{ag}));
            end
        end
    end
    
%!    % create structure array for "master gene list"
%!    genesnp = struct('geneID',uqids', ...
%!                     'snps',snps, ...
%!                     'rgindex',rgindex, ...
%!                     'chrn',num2cell(chrn), ...
%!                     'gene_start',num2cell(g_start), ...
%!                     'gene_end',num2cell(g_end), ...
%!                     'snp_start',num2cell(snp_start), ...
%!                     'snp_end',num2cell(snp_end) );
 
    % create index for compressing/ordering genes/snps
    geneHasSnps = find(cellfun(@(s) ~isempty(s),snps));
    snps = snps(geneHasSnps);
    firstsnps = cellfun(@(s) s(1),snps);
    [firstsnps,order] = sort(firstsnps);
    geneHasSnps = geneHasSnps(order);
    snps = snps(order);

    
    % create output struct by removing matrix fields
    matrix_fields = {'dat','orig','cbs','affy_calls'}; % TODO add more?
    if isfield(D,'matrix_fields')
        matrix_fields = unique([matrix_fields D.matrix_fields]);
    end
    G = rmfield_if_exists(D,matrix_fields);
    % remove some other fields we won't need
    G = rmfield_if_exists(G,{'marker','chr','pos','cM', ...
                            'score','cbs_rl','orig', ...
                            'history','origidx','gorigidx'});
    % collapse G rows from snp to gene size 
    G = reorder_D_rows(G,firstsnps);
    % add some gene fields
    G = add_D_field(G,'gene',{'geneID','snps','rgindex','chrn',...
                               'gene_start','gene_end'});
    G.geneID = uqids(geneHasSnps)';
    G.snps = snps;
    G.rgindex = rgindex(geneHasSnps);
    G.chrn = chrn(geneHasSnps);
    G.gene_start = g_start(geneHasSnps);
    G.gene_end = g_end(geneHasSnps);
    G.collapse_settings = options;
    
    %% collapse copy number data

    % loop across all copy number matrix fields in D
    for field = matrix_fields
        fld = char(field);
        if isfield(D,fld)
            fldmin = [fld 'min'];
            fldmax = [fld 'max'];
            G = add_D_field(G,'matrix',{fldmin,fldmax});
            G.(fldmin) = nan(length(firstsnps),size(D.(fld),2));
            G.(fldmax) = nan(length(firstsnps),size(D.(fld),2));
            % for speed, undo any SegArrays
            if isa(D.(fld),'SegArray')
                dat = full(D.(fld));
            end
            for i = 1:length(snps)
                G.(fldmin)(i,:) = min(dat(snps{i},:),[],1);
                G.(fldmax)(i,:) = max(dat(snps{i},:),[],1);
            end
        end
    end

