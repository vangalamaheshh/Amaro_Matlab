function [AC scramble_hp] = scrambleAchillesGenes( AC )
%scrambleAchillesGenes shuffle Achilles gene (row) data 
%   Randomly reorders hairpin data for genes among genes with the same 
%   number of hairpins in the Achilles data container AC. The hairpinIDs
%   are reordered along with their associated data so that the original
%   source of the data can be identified. The scramble_hp 2nd output identifies
%   the original location of the scrambled hairpins.

    % scramble the hairpin data by gene, 
    % permuting among genes w/equal hairpin representation
    [nhps ngenes] = size(AC.umap);
    scramble_gene = 1:ngenes;
    scramble_hp = 1:nhps;
    hp_counts = full(sum(AC.umap,1));
    % loop across sets of genes with same number of hairpins
    hps_per_gene = full(sum(AC.umap,1));
    for c=find(histc(hps_per_gene,1:max(hps_per_gene))>1)
        % reorder the genes in eaqch hairpin group
        gc = find(hp_counts == c);
        scramble_gene(gc) = gc(randperm(length(gc)));
    end
    % create reordered hairpin index
    for g = 1:ngenes
        ghps = find(AC.umap(:,g));
        scramble_hp(AC.umap(:,scramble_gene(g))) = ghps(randperm(length(ghps)));
    end 
    % scramble data
    AC.dat = AC.dat(scramble_hp,:);
    % scramble hairpin names to match data (for tracking results)
    AC.hairpinID = AC.hairpinID(scramble_hp);
end

