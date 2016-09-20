function Tc = condense_regions(T)
%
% condense_regions(T)
%
% Given region list T, a structure with required fields:
%    start
%    end
% and optional fields:
%    gene
%    chr
%    strand
% 
% Returns region list Tc
% in which all contiguous/overlapping/redundant regions have been consolidated.
% If gene, chr, and/or strand are supplied, will require that they match in order for consolidation to happen.
%
% Note: this function does not strip spliceform identifiers (.*)
% Use strip_spliceform_ids() beforehand for this purpose.
% Also it does not add flanking nucleotides; use
% flank_regions() for this. 
%
% Mike Lawrence 2008-05-15
%      modified 2012-04-03 to make optional these fields: gene, chr, strand


require_fields(T, {'start';'end'});
fake_fields = {};
if ~isfield(T,'gene'), T.gene = repmat({'notused'},slength(T),1); fake_fields{end+1} = 'gene'; end
if ~isfield(T,'chr'),  T.chr = repmat({'chrN'},slength(T),1); fake_fields{end+1} = 'chr'; end
if ~isfield(T,'strand'), T.strand = repmat({'+'},slength(T),1); fake_fields{end+1} = 'strand'; end

if isnumeric(T.chr), error('"chr" field should not be numeric'); end

[genes gi gj] = unique(T.gene);
ng = length(genes);

nt = slength(T);

Tc = [];
Tc.gene = cell(nt,1);
Tc.chr = cell(nt,1);
Tc.start = zeros(nt,1);
Tc.end = zeros(nt,1);
Tc.strand = zeros(nt,1);
nT = 0;

for gno=1:ng
     if ~mod(gno,100), fprintf('gene %d of %d\n', gno, ng); end

     gidx = find(gj==gno);
     chr = T.chr{gidx(1)};
     gene_min = min(T.start(gidx));
     gene_max = max(T.end(gidx));
     gene_len = gene_max - gene_min + 1;
     dna = zeros(gene_len,1);
     for t = 1:length(gidx)
       s = T.start(gidx(t));
       e = T.end(gidx(t));
       si = s - gene_min + 1;
       ei = e - gene_min + 1;
       dna([si:ei]) = T.strand{gidx(t)};    % fill in with "+" or "-"
     end

     pos = 1;
     strand = 0;
     while(1)
        % seek to start of next region
        while(pos<=gene_len)
          if dna(pos)
             strand = dna(pos);
             break
          end
          pos = pos + 1;
        end
        if pos > gene_len, break; end
        s = pos;
        % seek to end of this region
        while(pos<gene_len)
          if ~dna(pos+1), break; end
          if dna(pos+1) ~= strand
               error('gene %s abruptly switched strands!', genes{gno});
          end
          pos = pos + 1;
         end
        e = pos;
        nT = nT + 1;
        Tc.gene{nT} = genes{gno};
        Tc.chr{nT} = chr;
        Tc.strand(nT) = strand;
        Tc.start(nT) = s + gene_min - 1;
        Tc.end(nT) = e + gene_min - 1;
        pos = pos + 1;
     end
end

Tc = reorder_struct(Tc,1:nT);

Tc = rmfields(Tc,fake_fields);
