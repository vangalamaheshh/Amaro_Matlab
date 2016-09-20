function [ M ] = mini_driver(M, nperm, theta, output_file, ...
    report_file, image_dir, refseq_struct, first_idx, last_idx)
%The main function for the clustering program. Input parameters are:
%
%   Detailed explanation goes here

ensure_dir_exists(image_dir);
p_values = cell(slength(M{1}.gene), 5);
report = fopen(report_file, 'wt');
import('org.broadinstitute.cga.tools.seq.FixedWidthBinary');
forbidden = [-1,220,230,240];
nsets = length(M);
builds = cell(nsets, 1);

% PREPROCESS EACH DATASET

%if ~iscell(M), M = {M}; endw
%check_multiM_agreement(M);

setnames = cell(nsets,1);
for si=1:nsets
    setnames{si} = M{si}.name;
    % mark the mutations being used
    M{si}.mut.use_nonsilent = false(slength(M{si}.mut),1);
    M{si}.mut.use_nonsilent(M{si}.use_nonsilent) = true;
    % remove indel/null categories
    cidx = grepi('indel|null',M{si}.categ.name,1);
    M{si}.categ = reorder_struct_exclude(M{si}.categ,cidx);
    M{si}.mut = reorder_struct(M{si}.mut,~ismember(M{si}.mut.categ,cidx));
end



for si=1:nsets
    builds{si} = M{si}.build ;
end

hg19_idx = find(ismember(builds, 'hg19')==1, 1,'first');
hg18_idx = find(ismember(builds, 'hg18')==1, 1,'first');


%% Figure out what builds are present for the datasets used, 
% and prepare a bi indexing variable to loop over builds for each gene
% for the coordinate processing bellow

if isempty(hg18_idx)
  bi = 1; 
elseif isempty(hg19_idx) 
  bi = 2;
else 
  bi = [1, 2];
end
%keyboard


ghist = zeros(M{1}.ng,1);

for si=1:nsets
    ghist = ghist + histc(M{si}.mut.gene(M{si}.mut.use_nonsilent),1:M{1}.ng);
end


[tmpzz, gene_order] = sort(ghist,'descend');

if ~exist('first_idx','var'), first_idx = 1; end
if ~exist('last_idx','var'), last_idx = inf; end

for i = gene_order(5:length(gene_order))'  
    
    if i<first_idx || i>last_idx, continue; end
    
    %  keyboard
    name = M{1}.gene.name{i}; disp(name);
    
    % GENE characteristics
    starts = cell(2,1);
    ends = cell(2, 1);
    exons_matrix = cell(2, 1);
    genomic_array = cell(2,1);
    polyphen_track = cell(2,1);
    categories_track = cell(2,1);
    
    
    chr = M{1}.gene.chr(i);
    genelength = M{1}.gene.len(i);
    
    % CATEGORIES TRACKS
    
    fileCat_hg19 = FixedWidthBinary('/xchip/cga1/lawrence/db/hg19/context65/all.fwb');
    fileCat_hg18 = FixedWidthBinary('/xchip/cga1/lawrence/db/context65/all.fwb');
    % categories_track = populate_coverage(genomic_array{2}, fileCat, chr);
    
    
    % POLYPHEN TRACKS
    
    fileA_hg19 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg19/wigs/allA.cat.fwb');
    fileC_hg19 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg19/wigs/allC.cat.fwb');
    fileG_hg19 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg19/wigs/allG.cat.fwb');
    fileT_hg19 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg19/wigs/allT.cat.fwb');
    
    fileA_hg18 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg18/wigs/allA.cat.fwb');
    fileC_hg18 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg18/wigs/allC.cat.fwb');
    fileG_hg18 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg18/wigs/allG.cat.fwb');
    fileT_hg18 = FixedWidthBinary('/xchip/cga1/petar/reference_files/hg18/wigs/allT.cat.fwb');

    
    
    %% For each build, figure out the exons, the genomic coordinates as well as the polyphen and categories tracks %%
    %keyboard
    for bi = bi
        if bi == 1
            curr_idx = hg19_idx;
        else
            curr_idx = hg18_idx;
        end
        starts{bi} = M{curr_idx(1)}.cov.targ.start(M{curr_idx(1)}.cov.targ.gidx == i);
        ends{bi} = M{curr_idx(1)}.cov.targ.end(M{curr_idx(1)}.cov.targ.gidx == i);
        curr_starts = starts{bi};
        curr_ends = ends{bi};
        exons_matrix{bi} = [curr_starts'; curr_ends']';
        genomic_array{bi} = get_GenomicArray(exons_matrix{bi}, ...
            genelength);
        
        
        if bi == 1
            polyphen_track{bi} = populate_polyphen(genomic_array{bi}, fileA_hg19, fileC_hg19, fileG_hg19, fileT_hg19, chr);
            categories_track{bi} = populate_coverage(genomic_array{bi}, fileCat_hg19, chr);
        else
            polyphen_track{bi} = populate_polyphen(genomic_array{bi}, fileA_hg18, fileC_hg18, fileG_hg18, fileT_hg18, chr);
            categories_track{bi} = populate_coverage(genomic_array{bi}, fileCat_hg18, chr);
        end
        
        
    end
    
    
    
    %% Check the target list for errors (for sh-RNAs and such, mutliple target lists may come up as one)
    error_found = 0;
    for i = 1:length(exons_matrix)
      for j = bi
        if starts{bi}(1) > starts{bi}(length(starts{bi}))
          error_found = 1; 
          break 
        end
      end
    end

    if error_found
      fprintf('Error with target list retrieved for gene: %s\n', name);
      continue
    end
    % MUTATIONS in each tumor type
    
    mutations = cell(nsets,1);
    newbase = cell(nsets,1);
    num_mutations = zeros(nsets,1);
    mut_categories = cell(nsets,1);
    tot_num_mutations = 0;
    
    coverage_track = cell(nsets,1);
    maxCov = zeros(nsets,1);
    error_found = 0;
    keyboard
    for si = 1:nsets
        if ismember(builds{si}, 'hg19')
            curr_exons_matrix = exons_matrix{1};
            curr_genomic_array = genomic_array{1};
        else
            curr_exons_matrix = exons_matrix{2};
            curr_genomic_array = genomic_array{2};
        end
        midx = find(M{si}.mut.gene == i);
        midx = midx(M{si}.mut.use_nonsilent(midx) == true);
        midx = midx(strcmp(M{si}.mut.classification(midx), 'SNP') == 1);
        mutations{si} = M{si}.mut.start(midx);
     %   if any(mutations{si} <= min(starts{bi}
        newbase{si} = listmap(M{si}.mut.newbase(midx), {'A', 'C', 'G', 'T'});
        num_mutations(si) = size(mutations{si}, 1);
        tot_num_mutations = tot_num_mutations + num_mutations(si);
        %        keyboard
        mutations{si} = mutation_converter(curr_exons_matrix, genelength, mutations{si});
        mutations{si} = mutations{si}(mutations{si} ~= 0);
        mut_categories{si} = M{si}.mut.categ(midx);
        fileCov = FixedWidthBinary(M{si}.file.summed_cov_track);
        coverage_track{si} = populate_coverage(curr_genomic_array, fileCov, chr);
        maxCov(si) = max(coverage_track{si});
    end
    
    
    
    flag = 'both'; 
    if tot_num_mutations < 2
        flag = 'pph_only';
                
    else
        
        hg_19_bad = false; 
        hg_18_bad = false; 
        if ~(any(polyphen_track{1} ~= -1) & (sum(sum(polyphen_track{1} ~= forbidden(1) & polyphen_track{1} ~= forbidden(2) & polyphen_track{1} ~= forbidden(3) & polyphen_track{1} ~= forbidden(4)))  >= 0.55*genelength*4))
            fprintf('Problem with pph track\n');
            hg_19_bad = true; 
        end
        
        
        if ~(any(polyphen_track{2} ~= -1) & (sum(sum(polyphen_track{2} ~= forbidden(1) & polyphen_track{2} ~= forbidden(2) & polyphen_track{2} ~= forbidden(3) & polyphen_track{2} ~= forbidden(4)))  >= 0.55*genelength*4))
            fprintf('Problem with pph track\n');
            hg18_bad = true; 
        end
        
        if ~ismember(builds{si}, 'hg19')
            pph_to_use  = polyphen_track{2};
            cat_to_use = categories_track{2};
        else
            pph_to_use = polyphen_track{1};
            cat_to_use = categories_track{1};
        end
        
        
        if (~ismember(builds{si}, 'hg19') & hg18_bad) | (~ismember(builds{si}, 'hg18') & hg19_bad) & ~strcmp(flag, 'pph_only')
            flag = 'ks_only';
        end 
        
        if (~ismember(builds{si}, 'hg19') & hg18_bad) | (~ismember(builds{si}, 'hg18') & hg19_bad) & strcmp(flag, 'pph_only')
           p_ks = 1, p_pph = NaN, p_joint = 1; 
        end
        
        keyboard
        
        % GENERATE LIST OF THROWABLE MUTATIONS
        
        throwable = cell(nsets,1);
        mut_matrix = cell(nsets,1);
        for si=1:nsets
            Q = assign_65x4_to_categ_set(M{si}.categ);
            ncat = max(mut_categories{si});
            temp_cell_cont = cell(ncat,1);
            temp_cell_new = cell(ncat,1);
            temp_cell_total = cell(ncat,1);
            throwable{si} = cell(ncat, 1);
            for j = 1:ncat
                %keyboard
                [temp_cell_cont{j}, temp_cell_new{j}] = find(Q(:,j,:));
                temp_cell_total{j} = [temp_cell_cont{j}';temp_cell_new{j}']';
                throwable{si}{j} = populate_cont_pairs(coverage_track{si}, cat_to_use, temp_cell_total{j},pph_to_use);
            end
            if ~isempty(mutations{si})
                mut_matrix{si} = [mutations{si} cat_to_use(mutations{si}) coverage_track{si}(mutations{si}) pph_to_use(mutations{si}, 1)];
            end
        end
        %        keyboard
        
        keyboard
        
        
        % PERFORM PERMUTATIONS
        
        [ks_stat_array, p_ks, p_pph, p_joint, ci_ratio] = generate(genelength, nsets, setnames, tot_num_mutations, nperm, theta, mut_matrix, pph_to_use, mut_categories, ...
            mutations, chr, name, newbase, coverage_track, throwable, maxCov, image_dir, flag);
        
    end
    
    
    
    if p_ks ~= 1
        %fprintf(report, 'name = %s\n', name);
        fprintf(report, 'i = %f\n', i);
        fprintf(report, 'mutations = \n');
        all_mutations = cat(1,mutations{:});
        for j = 1:size(all_mutations,1)
            fprintf(report, '%d\n', all_mutations(j));
        end
        fprintf(report, 'len = %d\n', genelength);
        fprintf(report, 'p_kst = %f\n', p_ks);
        fprintf(report, 'p_pph = %f\n', p_pph);
        fprintf(report, 'p_joint = %f\n', p_joint);
        fprintf(report, '\n');
    end
    
    %    keyboard
    
    p_values{i,2} = p_ks;
    p_values{i,3} = p_pph;
    p_values{i,4} = p_joint;
    p_values{i,5} = ci_ratio;
    %    keyboard
    %  catch
    %    fprintf('Error with gene %s\n', name);
    %  end
end
fclose(report);

tmp = nan(size(p_values,1),1);
for i=1:length(tmp), if ~isempty(p_values{i,4}), tmp(i) = p_values{i,4}; end, end
[tmpzz,ord] = sort(tmp);
p_values = p_values(ord,:);


result_file = output_file;
save(result_file, 'p_values');

