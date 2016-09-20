function genelist_survey_FDgrid(inname,outname,build)
%
% Given a list
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3N_breakdown.m
%
% for each target transcript,
% determines if the transcript is a "perfect" ORF,
% and if so, computes Nsil + Nmis + Nnon = 3N
% for each base category, and outputs table.
%
% Mike Lawrence 2008/04/29
%
% modified 2008-05-21
%
% 1. to look at taking into account weights of mutational contexts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

addpath ~/CancerGenomeAnalysis/trunk/matlab/
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
addpath ~/CancerGenomeAnalysis/trunk/matlab/seq
addpath ~/CancerGenomeAnalysis/trunk/matlab/mike

project = 'TCGA';

if strcmp(project, 'TCGA')
  build = 'hg18';
  inname = '/xchip/tcga/gbm/analysis/lawrence/cov/all_targets.fixed_names.txt';
  outname = '/xchip/tcga/gbm/analysis/lawrence/cov/TCGA_breakdown.txt';
  ratio_outname = '/xchip/tcga/gbm/analysis/lawrence/cov/TCGA_ratio.txt';
elseif strcmp(project, 'TSP')
  build = 'hg17';
  inname = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_targets.txt';
  outname = '/xchip/tcga/gbm/analysis/lawrence/tsp/TSP_breakdown.txt';
  ratio_outname = '/xchip/tcga/gbm/analysis/lawrence/cov/TSP_ratio.txt';
  TSP_wash250_only = true;
end

output_ratio_file = true;

% WEIGHTS

%     1. A
%     2. T
%     3. C in CpG
%     4. C in TpC but not TpCpG
%     5. other C
%     6. G in CpG
%     7. G in GpA but not CpGpA
%     8. other G

weights_TSP = [0.36 0.36 2.84 1.4 1.4 2.84 1.0 1.0];    % from TSP
weights_TCGA = [0.45 0.45 3.2 0.83 0.83 3.2 1.2 1.2];   % from TCGA
weights_none = [1 1 1 1 1 1 1 1];                       % unweighted 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read uncondensed target list
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Reading target list\n');

if strcmp(project,'TCGA')
  table = read_table(inname, '%s%s%s%f%f%s%s', char(9), 1, 'Whitespace', '\b\r');
  target.gene = table.dat{1};
  target.center = table.dat{2};
  target.chr = table.dat{3};
  target.start = table.dat{4};
  target.end = table.dat{5};
  target.strand = table.dat{6};
  target.phase = table.dat{7};
  clear table;

  % Note: for this analysis, we will NOT collapse alternate spliceforms!

elseif strcmp(project, 'TSP')
  target = tab2struct(inname, '%s%s%s%s%s%f%f%s%s%s%s%s%s%s%s');
  target = rename_field(target, 'stop', 'end');
  target = rename_field(target, 'ref', 'chr');
  target = rename_field(target, 'name', 'gene');

  % try collapsing alternate spliceforms
  target.gene = strip_extensions(target.gene);

  wash250 = tab2struct('/xchip/tcga/gbm/analysis/lawrence/tsp/was250.txt');
  aliases = tab2struct('/xchip/tcga/gbm/analysis/lawrence/tsp/aliases.txt');
  wash250.name = apply_aliases(wash250.name, aliases, 'ignore_extensions');
  target.gene = apply_aliases(target.gene, aliases, 'ignore_extensions');
end

ntarg = length(target.gene);

% convert chromsome identifiers
target.chr = convert_chr(target.chr);

% convert genes to gene numbers
[gene_name gi target.g] = unique(target.gene);
ng = length(gene_name);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  COMPUTE 3N BREAKDOWN STATISTICS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

silent=1;

glen = zeros(ng,1);
stops = zeros(ng,6);
perfect_frame = zeros(ng,1);

sil_tot = zeros(8,1);
non_tot = zeros(8,1);
mis_tot = zeros(8,1);

out = fopen(outname, 'wt');

fprintf(out, 'gene');
for i=1:8, fprintf(out, '\tsil%d',i); end
for i=1:8, fprintf(out, '\tmis%d',i); end
for i=1:8, fprintf(out, '\tnon%d',i); end
fprintf(out,'\n');

if output_ratio_file
  ratio_out = fopen(ratio_outname,'wt');
  fprintf(ratio_out, 'gene\tlength\tunweighted\tTCGA_weights\tTSP_weights\n');
end

total = 0;

for g=1:ng          % for each gene

    gname = gene_name{g};
 
    if strcmp(project, 'TCGA')
      if strncmp(gname, 'ConsReg', 7), continue; end
      if strncmp(gname, 'hsa-', 4), continue; end

      %    if use_condensed_target_list
      %       if strncmp(gname, 'DAZ', 3), continue; end   % skip weirdo Y-gene DAZ
      %       if strncmp(gname, 'SPANXA', 6), continue; end  % also has strand-switching
      %    end
    elseif strcmp(project, 'TSP')
      if TSP_wash250_only
        g2 = regexprep(gname, '\..*', '');
        if ~ismember(g2,wash250.name), continue; end
      end
    end

    gene_targets = find(target.g == g);

    gene_chr = target.chr(gene_targets(1));
    if isnan(gene_chr) continue; end         % skip genes on "weird" chromosomes

    fprintf('gene %d/%d %s\n', g, ng, gname);

    % sort by start position and remove redundant exons
    [u ui uj] = unique(target.start(gene_targets));
    gene_targets = gene_targets(ui);

    % if gene is on (-) strand, reverse the order of the targets

    strand = target.strand(gene_targets(1));
    if strcmp(strand, '-'), gene_targets = flipud(gene_targets); end

    gene_dna = [];
    gene_context = [];
 
    for tno = 1:length(gene_targets)     % Note: a gene's targets are
                                         % guaranteed to be non-overlapping,
                                         % after the analysis in condense_targets.m
      t = gene_targets(tno);

      if ~strcmp(strand,target.strand(t)), error('Gene switches strands!\n'); end

      ts = target.start(t);
      te = target.end(t);
      tlen = te - ts + 1;
      glen(g) = glen(g) + tlen;

      seq = genome_region(gene_chr, ts-1, te+1, build);   % retrieve 1-nt flank for computing context
      if strcmp(strand, '-'), seq = my_seqrcomplement(seq); end

      context = survey_contexts(seq);
      context = context(2:end-1);
      seq = seq(2:end-1);        % remove 1-nt flank

      gene_dna = [gene_dna seq];
      gene_context = [gene_context context];

    end   % next target in this gene

    % count stops in each frame

    for f=1:6
        stops(g,f) = count_stops(gene_dna,f);
    end        

    % finished processing this gene:

    ATGstart = strcmpi(gene_dna(1:3),'atg');
    TERMend =  strcmpi(gene_dna(end-2:end),'tag') | strcmpi(gene_dna(end-2:end),'tga') ...
      | strcmpi(gene_dna(end-2:end),'taa');

    lenmod3 = mod(glen(g),3);

% old condition (very stringent):
%    is_perfect(g) = (lenmod3==0)&(stops(g,1)==1)&(ATGstart==1)&(TERMend==1);

% new condition (much less stringent:
%    if sum(stops(g,1:3)<=1)>1, error('More than one perfect frame!\n'); end
    for f=1:3
      if stops(g,f)<=1
        perfect_frame(g) = f;
        break;
      end
    end

    if perfect_frame(g)
       % examine each possible mutation in silico to see what kind of mutation it causes
       M = survey_all_mutations(gene_dna(perfect_frame(g):end));
       sil = zeros(8,1); mis = zeros(8,1); non = zeros(8,1);
       for c=1:8
         pos = find(gene_context==c);
         pos(pos>size(M,2)) = [];    % in case of lenmod3~=0
         sil(c) = sum(M(1,pos));
         mis(c) = sum(M(2,pos));
         non(c) = sum(M(3,pos));
       end
    else

       % nonperfect
       sil=-ones(8,1);
       mis=-ones(8,1);
       non=-ones(8,1);
    end

    % compute R, expected of silent to nonsilent mutations

    if perfect_frame(g)
      R_unweighted = (weights_none * sil) / (weights_none * (mis+non));
      R_TCGA = (weights_TCGA * sil) / (weights_TCGA * (mis+non));
      R_TSP = (weights_TSP * sil) / (weights_TSP * (mis+non));

      sil_tot = sil_tot + sil;
      mis_tot = mis_tot + mis;
      non_tot = non_tot + non;
    else
      R_unweighted = -1;
      R_TCGA = -1;
      R_TSP = -1;
    end

    fprintf(out, '%s', gene_name{g});
    for i=1:8, fprintf(out, '\t%d', sil(i)); end
    for i=1:8, fprintf(out, '\t%d', mis(i)); end
    for i=1:8, fprintf(out, '\t%d', non(i)); end
    fprintf(out,'\n');

    if output_ratio_file
      fprintf(ratio_out, '%s\t%d\t%f\t%f\t%f\n', gene_name{g}, glen(g), R_unweighted, R_TCGA, R_TSP);
    end

    total = total + 1;

end   % next gene

fclose(out);

if output_ratio_file
  fclose(ratio_out);
end

sil_tot
mis_tot
non_tot

fprintf('Total perfect transcripts = %d / %d\n', sum(perfect_frame>0), total);

