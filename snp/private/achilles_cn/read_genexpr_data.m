function X = read_genexpr_data(fname)
%READ_GENEXPR_DATA read in Achilles-formatted data table.
%
%   X = READ_GENEXPR_DATA(FNAME) reads an Achilles data table into
%   the structure X. The input data are in a tab-delimited text file
%   with a header. The first column is of gene IDs, the next
%   column is of probe IDs, the remaining columns are z-scored
%   essentiality data. 
%
%   The return value contains the achilles dataset:
%       sampleID - contains sample identifiers (cell-string)
%       hairpinID - contains the unique probe ID
%       geneID - names the gene the hairpin was meant to knock down
%       dat - hairpinID x sampleID array of Z-scores (doubles) 
%       ugene - unique subset of geneID
%       umap - logical map of hairpins (row) to genes (column)


    % read the header to get sample IDs (every other field starting at 3)
    [hdr,fid] = read_dlm_file(fname,char(9),3);
    X = struct;
    X.sampleID = hdr{1}(3:2:end);
    % read the data
    Nsamples = length(X.sampleID);
    form = ['%s %s' repmat(' %f %*s',1,Nsamples)];
    datacells = textscan(fid,form,'delimiter',char(9));
    % first column gene ID, second column probe ID
    X.geneID = datacells{1};
    X.hairpinID = datacells{2};
    X.dat = [datacells{3:end}];

    % generate list of unique genes and mapping of "hairpins" to it
    X.ugene = unique(X.geneID);
    [~,gx,ux] =match_string_sets_hash(X.geneID,X.ugene);
    X.umap = sparse(gx,ux,true);
%!    X.umap = logical(match_string_sets_hash(X.geneID,X.ugene));
