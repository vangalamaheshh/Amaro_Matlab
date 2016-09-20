function V = load_vog_muts()
% load all Vogelstein mutations

dirname = '/xchip/tcga/gbm/analysis/lawrence/vog/';

% breast and colorectal

fprintf('Loading breast and colorectal data\n');
brco = load_struct([dirname 'brco/tableS3_liftover.txt']);
brco = rmfield(brco,{'sift_score','log_RE_value','LS_SNP_score'});

% GBM

fprintf('Loading GBM data\n');
gbm1 = load_struct([dirname 'gbm/tableS3_liftover.txt']);
% remove hypermutated tumor
gbm1 = reorder_struct(gbm1,~strcmpi(gbm1.tumor,'Br27P'));
gbm1 = rmfield(gbm1,{'hg18_site','LS_MUT_score'});
gbm1.tumor_type = repmat({'GBM'},slength(gbm1),1);
gbm1.screen = repmat({'Discovery'},slength(gbm1),1);
gbm2 = load_struct([dirname 'gbm/tableS4_liftover.txt']);
gbm2.tumor_type = repmat({'GBM'},slength(gbm2),1);
gbm2.screen = repmat({'Prevalence'},slength(gbm2),1);

% pancreatic

fprintf('Loading pancreatic data\n');
panc1 = load_struct([dirname 'panc/tableS3_liftover.txt']);
panc1.tumor_type = repmat({'Pancreatic'},slength(panc1),1);
panc1.screen = repmat({'Discovery'},slength(panc1),1);
panc1 = rmfield(panc1,'LS_MUT_score');
panc2 = load_struct([dirname 'panc/tableS4_liftover.txt']);
panc2.tumor_type = repmat({'Pancreatic'},slength(panc2),1);
panc2.screen = repmat({'Prevalence'},slength(panc2),1);

% combine all data

fprintf('Combining all data\n');
V.mut = combine_structs({brco,gbm1,gbm2,panc1,panc2});
[V.gene.name tmp V.mut.gene] = unique(V.mut.gene);
V.mut.site = V.mut.hg18;
V.ng = length(V.gene.name);
V.mut.start = str2double(V.mut.start);
V.mut.end = str2double(V.mut.end);
V.mut.chr = convert_chr(V.mut.chr);
