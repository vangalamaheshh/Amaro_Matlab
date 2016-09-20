function [rg version] = make_gaf_refgene(gaf_path)
% MAKE_GAF_REFGENE create referenc genome from GAF file
%
%   [RG VERSION] = make_gaf_refgene(GAF_PATH)
%
% Returns a reference genome struct in RG and a version string of
% make_gaf_refgene in VERSION. GAF_PATH is a path to the GAF text file,
% currently named something like TCGA.hgXX.<monthyear>.gaf

% version of GAF maker
version = '1.0.1';

% select the FeatureType values to keep
keep_types = {'gene','pre-miRNA'};
%!keep_types = {'gene','miRNA','pre-miRNA'};
%! NOTE: the miRNA types are problematic, some have no genomic position,
%! some have conflicting positions under CompositeCoordinates and GeneLocus 

%% open and read the GAF file
fid = fopen(gaf_path);
header = fgetl(fid);
hdr_fields = regexp(header,char(9),'split');
ncols = length(hdr_fields);
fprintf('%d columns detected in ''%s'', loading data...\n',ncols,gaf_path);
row_format = repmat('%s',1,ncols);
gaftxt = textscan(fid,row_format,'Delimiter',char(9),'Bufsize',2e5);
fclose(fid);
nrows = length(gaftxt{1});
fprintf('loaded %d rows of data\n\n',nrows);

%% check out how many of what are in the file
% third column is feature type
assert(strcmp('FeatureType',hdr_fields{3}));
fprintf('examining FeatureType content of file...\n');
feat_types = unique(gaftxt{3});
% also flag the rows we want to keep
for k = keep_types
    assert(ismember(k,feat_types));
end
keepers = false(nrows,1);

for i = 1:length(feat_types)
    feat_type = feat_types{i};
    featurers = strcmp(feat_type,gaftxt{3});
    fprintf('%d\t%s',sum(featurers),feat_type);
    if ismember(feat_type,keep_types)
        keepers = keepers | featurers;
        fprintf(' **');
    end
    fprintf('\n');
end

%% reduce the data to the features we want
fprintf('\n** reducing features from %d to %d of interest\n',length(keepers),sum(keepers));
for i=1:length(gaftxt)
    gaftxt{i} = gaftxt{i}(keepers);
end

%% get gene location from GenLocus
assert(strcmp('GeneLocus',hdr_fields{17}));
split_locus = regexp(gaftxt{17},'^(chr.[^:]+):([0-9]+)-([0-9]+):([+-])$','tokens');
% NOTE: regexp pattern chosen to exclude multiple loci separated by semicolons

%% get gene location from CompositeCoordinates
assert(strcmp('CompositeCoordinates',hdr_fields{15}));
split_coord = regexp(gaftxt{15},'^(chr.+):([0-9]+)[0-9,-]*-([0-9]+):([+-])$','tokens');
% NOTE: regexp designed to collapse a disjoint range into start and end location  

%% identify data with good positions
pos_from_locus = 0<cellfun(@length,split_locus);
pos_from_coord = 0<cellfun(@length,split_coord);
dual_pos = pos_from_locus & pos_from_coord;
no_pos = ~(pos_from_locus|pos_from_coord);
fprintf('Detected %d feature(s) with no genomic position\n',sum(no_pos));

% find all locations where positions are inconsistent
pos_conflicts = false(size(split_locus));
for i=1:length(split_locus)
    pos_conflicts(i) = ~isequal(split_locus{i},split_coord{i});
end
pos_conflicts = pos_conflicts & dual_pos;
fprintf('Detected %d genomic position conflict(s)\n',sum(pos_conflicts));

%% keep only features with genomic positions that are consistent
keepers = ~(no_pos|pos_conflicts);
fprintf('Retaining %d features (genes) with good genomic positions\n',sum(keepers));
% shrink table
gaftxt = horzcat(gaftxt{:});
gaftxt = gaftxt(keepers,:);
%combine parsed locations
split_locus(~pos_from_locus) = split_coord(~pos_from_locus);
split_locus = split_locus(keepers);

%% split feature ID into genes and identifiers
assert(strcmp('FeatureID',hdr_fields{2}));
split_featid = regexp(gaftxt(:,2),'\|','split');
name = cellfun(@(x) x{1},split_featid,'Uniform',false);
geneid =  cellfun(@(x) x{2},split_featid,'Uniform',false);
ismulti = 3==cellfun(@length,split_featid);
genstance = repmat({''},size(split_featid));
genstance(ismulti) = cellfun(@(x) x{3},split_featid(ismulti),'Uniform',false);

%% create minimum unique symb
symb = name;
% differentiate genes with same name by adding |id
[~,m1,m2]=match_string_sets_hash(symb,symb);
dups = full(sum(sparse(m1,m2,true))>1);
symb(dups) = strcat(symb(dups),'|',geneid(dups));
% further differentiate duplicate name|id s by adding kofn info
[~,m1,m2]=match_string_sets_hash(symb,symb);
dups = full(sum(sparse(m1,m2,true))>1);
symb(dups) = strcat(symb(dups),'|',genstance(dups));

%% create locus_id
% numeric value for compatibility purposes
%!!! TODO stop using locus_id
lex_geneid = regexp(geneid,'([^0-9]*)([0-9]+)','tokens');
lex_geneid = vertcat(lex_geneid{:});
prefixed = cellfun(@(x) ~isempty(x{1}),lex_geneid);
locus_id = cellfun(@(x) str2double(x{2}),lex_geneid);
locus_id(prefixed) = -locus_id(prefixed);

%% genomic coordinates
split_locus = vertcat(split_locus{:});
chr = cellfun(@(x) x{1},split_locus,'Uniform',false);
start = cellfun(@(x) str2double(x{2}),split_locus);
stop = cellfun(@(x) str2double(x{3}),split_locus);
strand = cellfun(@(x) strcmp('+',x{4}),split_locus);
chrn = chromosome2num(chr);

%% package it into a struct
% NOTE: some fields are cast to int32to be compatible with GISTIC
rg = struct('symb',symb,...
            'gene',name,...
            'locus_id',num2cell(int32(locus_id)),...
            'gaf_name',gaftxt(:,2),...
            'chr',chr,...
            'strand',num2cell(strand),...
            'start',num2cell(int32(start)),...
            'end',num2cell(int32(stop)),...
            'type',gaftxt(:,3),...
            'chrn',num2cell(chrn));
%!!! TODO 'gene' field replaced by descriptive name
%!!! TODO supply 'cds_start' and 'cds_end' fields