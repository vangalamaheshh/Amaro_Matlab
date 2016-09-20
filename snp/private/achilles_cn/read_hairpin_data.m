function AC = read_hairpin_data(fname)
%READ_HAIRPIN_DATA read in Achilles-formatted data table.
%
%   AC = READ_HAIRPIN_DATA(FNAME) reads an Achilles data table into
%   the structure AC. The input data are in a tab-delimited text file
%   with a header. The first column is of hairpin IDs, the next
%   column is of gene IDs, the remaining columns are z-scored
%   essentiality data. 
%
%   The return value contains the achilles dataset:
%       sampleID - contains sample identifiers (cell-string)
%       hairpinID - contains the unique hairpin ID
%       geneID - names the gene the hairpin was meant to knock down
%       dat - hairpinID x sampleID array of Z-scores (doubles) 
%       ugene - unique subset of geneID
%       umap - logical map of hairpins (row) to genes (column)

    % Use read_mit_gct_file() to read the file, but rename the
    % fields 'gacc' and 'gdesc' to 'hairpinID' and 'geneID' to 
    % make sense for Achilles 
    AC = read_mit_gct_file(fname);
    AC.hairpinID = AC.gacc;
    AC = rmfield(AC,'gacc');
    AC.geneID = AC.gdesc;
    AC = rmfield(AC,'gdesc');

    % make cell array sampleID field
    AC.sampleID = cellstr(AC.sdesc);

    % generate list of unique genes and mapping of hairpins to it
    AC.ugene = unique(AC.geneID);
    AC.umap = logical(match_string_sets_hash(AC.geneID,AC.ugene)); %! TODO fix this
