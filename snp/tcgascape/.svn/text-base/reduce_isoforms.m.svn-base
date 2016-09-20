function redgene = reduce_isoforms(rg, D)
%MAKE_REDGENE reduce genes to uniquer names largest isoform

% derived from first part of reduce_to_genes - could potentially be
% factored out

% compress genes with same identifier in reference genome
ids = {rg.symb};
uqids = unique(ids);
[~,m1,m2] = match_string_sets_hash(ids,uqids); %!
gmatch = sparse(m1,m2,true); %!
%!  [gmatch,m1,m2]=match_string_sets_hash(ids,uqids);
% m1 indexes ids
% m2 corresponds ids, indexes uqids
nred = length(uqids);

% full genome index selecting nothing
no_genome = SegArray.fromSegments(1,size(D.dat,1),1,false);

redgene = struct(...
    'symb',     uqids, ...
    'locus_id', Inf, ...
    'chrn',     0, ...
    'start',    Inf, ...
    'end',      0, ...
    'snps',     repmat({[]},1,nred), ...
    'segs',     repmat(no_genome,1,nred), ...
    'iso_idx',  repmat({[]},1,nred), ...
    'iso_segs', repmat({[]},1,nred) ...
);

% loop over unique gene symbols
for ag = 1:nred
    if mod(ag,1000)==0
        verbose('reduce_isoforms: processing gene %d of %d',30,ag,length(uqids));
    end
    matches = find(gmatch(:,ag));
    if ~isempty(matches)
        g_start = Inf;
        g_end = 0;
        g_chrn = 0;
        g_locus = Inf;
        g_snps = [];
        iso_segs = repmat({},1,length(matches));

        % loop over all matches to unique gene ID
        keepmatch = true(1,length(matches));
        for j = 1:length(matches)
            gene = rg(matches(j));
%!      for k = matches'
            if g_chrn == 0
                g_chrn = gene.chrn;
            elseif g_chrn ~= gene.chrn;
                verbose('Conflicting chromosomes (%d and %d) for gene %s !',30, ...
                        g_chrn,gene.chrn,gene.symb);
                if gene.chrn < g_chrn
                    g_chrn = gene.chrn;
                    % discard all previous matches
                    keepmatch = j < (1:length(matches));
                    g_snps = [];
                    g_start = Inf;
                    g_end = 0;
                    g_locus = Inf;
                else
                    % discard this match
                    keepmatch(j) = false;
                    continue;
                end
            end
            g_start = min(g_start,gene.start);
            g_end = max(g_end,gene.end);
            if g_end > g_start
                % combine markers for the same gene
%!              if isfield(rg,'snps')
                g_snps = union(g_snps,gene.snps);
%!              else
%!                  g_snps = union(g_snps,find_snps(D,g_chrn,g_start,g_end,1));
%!              end
            else
                verbose('Bad gene data for %s !',30,gene.symb);
            end
            g_locus = min(g_locus, gene.locus_id);
        end % loop over isoforms
        
        matches = matches(keepmatch);
        redgene(ag).iso_idx = matches;
        
        if ~all(diff(g_snps)==1)
            %!!! remove all but one contigous segment
            keepsnps = g_snps(1:find(diff(g_snps)>1,1,'first'));
            keepmatch = keepmatch & cellfun(@(x) ~isempty(intersect(keepsnps,x)), {rg(matches).snps});
            g_start = min([rg(keepmatch).start]);
            g_end = max([rg(keepmatch).end]);
            g_locus = min([rg(keepmatch).locus_id]);
        end
        
        % build logical index for isoform markers among all gene markers 
        for i = 1:length(matches)
            iso_segs{i} = ismember(g_snps,rg(matches(i)).snps);
        end
       
        % create segarray to cover the isoforms
        gaps = find(1~=diff(g_snps));
        B = [min(g_snps);g_snps(gaps+1)];
        E = [g_snps(gaps);max(g_snps)];
        redgene(ag).segs = SegArray.fromSegments(B,E);      
        
        % move local data into array
        redgene(ag).start = g_start;
        redgene(ag).end = g_end;
        redgene(ag).chrn = g_chrn;
        redgene(ag).locus_id = g_locus;
        redgene(ag).snps = g_snps;
        redgene(ag).iso_segs = iso_segs;
    end % if matches
end % loop over unique gene symbols

% remove genes without any markers
geneHasSnps = cellfun(@(s) ~isempty(s),{redgene.snps});
redgene = redgene(geneHasSnps);
