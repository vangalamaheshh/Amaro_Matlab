function [G SI CL Qs] = process_cle_data(input,cache,options,rg,cyto)

% TODO head documentation

% script for processing segmentsd CLE data into cached CN data and
% segment m-files schum 27 May 2010

% the outputs are essentially CL, Qs, and cyto
% the inputs are numerous

%% get cleaned copy number data
options = impose_default_value(options,'join_segment_size',10);
options = impose_default_value(options,'use_segarray',true);
options = impose_default_value(options,'remove_X',true);
options = impose_default_value(options,'remove_nans',true);
options = impose_default_value(options,'gene_minmax',false);
[CL SI] = get_clean_D_from_seg(cache.D_file,input,options);

% add cytoband information
CL = add_cyto(CL,cyto);

%% get ziggurat analysis of copy number data
if 0 %%%%%!!!!!
options.cn_cap = 1.5;
options.D_file = cache.D_file;
[Qs, CL] = get_ziggurat_segments(cache.Qs_file,CL,cyto,options);

% convert copy number data back to log units
% (ziggurat analysis converts them to linear)
if ~CL.islog
    CL.dat = log2(CL.dat+2)-1;
    CL.islog = true;
end

%% reconstruct broad and focal genomes

% (The segments of each cell element of Qs are summed into
% copy number array in a corresponding cell of the returned
% genome. Like Qs, the genome data units are copy
% number-2. Deletions are represented as positive numbers in this
% schema, e.g. total deletion from 2 copies has value +2.) 

% reconstruct the broad genome
verbose('=> Reconstructing broad genome.',10);
options.broad_or_focal = 'broad';
broad_genome = reconstruct_genomes(Qs,options);
% add field for broad events
CL = add_D_field(CL,'matrix','bat');
CL.bat = log2(2+broad_genome.amp-broad_genome.del) - 1;
clear broad_genome;

% reconstruct the focal genome
verbose('=> Reconstructing focal genome.',10);
options.broad_or_focal = 'focal';
focal_genome = reconstruct_genomes(Qs,options);
% add field for focal events
CL = add_D_field(CL,'matrix','fat');
CL.fat = log2(2+focal_genome.amp-focal_genome.del) - 1;
clear focal_genome;
end %%%%%!!!!
%% collapse SNP-level data to gene level data
if exist(cache.G_file,'file')
    G = load_D(cache.G_file);
else
    verbose('=> Constructing gene-level copy number data from SNP-level data.',10);
    if options.gene_minmax
        G = snpD_to_minmax_geneG(CL,rg,'symb');
        save_D(cache.G_file,G,'-v7.3');
    else
        G = snpD_to_geneG(CL,rg,'symb');
        save_D(cache.G_file,G,'-v7.3');
    end
end
