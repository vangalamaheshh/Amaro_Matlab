% make the gene summary tables which were deliberately removed from the
% script for the 4/10 and 4/16 runs
% hacked from tumorscape web file generator
% 2012-04-17 Steve Schumacher DFCI

clear
dbstop if error;
set_verbose_level(30);

%% parameters
% minimum number of samples required for inclusion in tumorscape
MIN_SAMPLES = 40;
%!!! need external reference genome for cyto - get from gistic_inputs
rg_file = '/xchip/gistic/variables/hg19/hg19_with_miR_20120227.mat';
%!rg_file = '/xchip/gistic/variables/hg18/hg18_with_miR_20080407.mat';

%! these parameters should come from gistic_inputs.mat from a GISTIC output directory
thresh = [0.3 0.3];
high_thresh = [1.0 1.0];
broad_len_cutoff = 0.5;

use_segarray = true;

%!addpath ('~/CancerGenomeAnalysis/trunk/matlab/snp/tcgascape/');

% date file extension for lookup tables
datext = '20120416';

%% working directories
tcgascape_dir = '/xchip/gistic/tcgascape/';    % all runs
run_dir = [tcgascape_dir 'tcgascape_120416/']; % this run
web_dir = [run_dir 'webfiles/'];  % flat files
subtype_dir = [run_dir 'gistic_results/'];     % inputs

tic;

%% create output directories
if ~exist(run_dir, 'dir')
    mkdir(run_dir);
end
if ~exist(web_dir, 'dir')
    mkdir(web_dir);
end

save_analysis_file = [run_dir 'webfile_data.' datext '.mat'];

% load cytoband information
load(rg_file)
clear rg

% load the master data, extract amplification and deletion events
D = load_D([subtype_dir 'all_cancers/D.cap2.all_cancers.mat']);
if ~use_segarray && strcmp(class(D.dat),'SegArray');
    verbose('Converting D.dat from SegArray to full array',20);
    D.dat = full(D.dat);
end

% TODO nail down what this file needs to be
QA = D.Qs.amp;
QD = D.Qs.del;

%% construct list of tumor types

root_type = 'all_cancers';

% load cancer type hierarchy and name mappings
[cancer_namemap cancer_treegen] = tcga_cancer_types();

% find all the unique cancer types in the dataset
sample_types = {D.sis.gcmtype};
typecounts = count_cancer_types(root_type,sample_types,cancer_treegen);
unique_types = keys(typecounts);

% output the counts
counts = cell2mat(values(typecounts,unique_types));
verbose('Tumorscape run - %s',10,datext);
for i = 1:length(unique_types)
    verbose('%s\t%d',10,unique_types{i},counts(i));
end

% remove cancer types without enough samples
keep = counts >= MIN_SAMPLES;
unique_types = unique_types(keep);
ntypes = length(unique_types);
count = counts(keep);
root_type_idx = find(strcmp(root_type,unique_types),1);
non_rootypes = setdiff(1:ntypes,root_type_idx);

% create a cell array of sample indexes for each type
al_idx = values(index_cancer_types(root_type,sample_types,cancer_treegen),unique_types);

% determine base types or "independent subsets"
isleaf = leaf_cancer_types(root_type,cancer_treegen);
isubsets = cell2mat(values(isleaf,unique_types));

%% load cancer type specific gistic results

sub_regs = cell(1,ntypes);    % wide peak regions for each subtype
sub_q = cell(1,ntypes);       % q value for each subtype
gene_regs = cell(1,ntypes);   % genes for each subtype region


for j=1:ntypes
  cur_type = char(unique_types(j));
  verbose('loading results for ''%s''',20,cur_type);
  % information about peaks
  load([subtype_dir cur_type '/wide_peak_regs.' cur_type '.mat']);
  sub_regs{j} = wregs{1}; 
  %! (wide_peak struct now stored multiply in cell array)
  % gene gistic statistical results
  load([subtype_dir cur_type '/gene_stats.' cur_type '.mat'])
  % make sure reference genomes are consistent across disease types
  if j > 1
      assert(length(gene_gistic_rg) == length(gene_regs{1}));
      assert(all(strcmp({gene_gistic_rg.symb},{gene_regs{1}.symb})));
  end
  gene_regs{j} = gene_gistic_rg;
  % snp level q values
  load([subtype_dir cur_type '/orig_stats.' cur_type '.mat'])
  sub_q{j} = q;
end

%% translate TCGA names to user-visible names
% (they were used to reference subdirectory names, new names will be those
% in tumorscape database)
user_types = values(cancer_namemap,unique_types);
newtypes = values(cancer_namemap,{D.sis.gcmtype});
for k=1:length(newtypes)
  D.sis(k).gcmtype = newtypes{k};
end

%% reduce genes to the largest isoform
urg = reduce_isoforms(gene_regs{1},D);

%{
%% Keep provisional genes - Remove random genes
goodchrn = find([urg.chrn] <= 23);
urg = urg(goodchrn);
[unique_symbs ui uj] = unique({urg.symb}); 
%}


%% Fix issues where there are more than one deletion peak with exactly same boundaries

verbose('fixing duplicate peaks',20);
for j=1:length(sub_regs)
  disp(j)
  cur_regs = sub_regs{j};
  cur_reg = cur_regs{2};
  rm = [];
  for i=1:length(cur_reg)
    conflict = find([cur_reg.peak_wide_st] == cur_reg(i).peak_wide_st & ...
           [cur_reg.peak_wide_en] == cur_reg(i).peak_wide_en);
    if any(conflict ~= i)
      [m mi] = min([cur_reg(conflict).resid_qv]);
      rm = unique([rm conflict(setdiff(1:length(conflict),mi))]);
    end
  end
  cur_reg = cur_reg(setdiff(1:length(cur_reg),rm));
  cur_regs{2} = cur_reg;
  sub_regs{j} = cur_regs;
end

% add gene lists (exclude partial hits for deletions)
for j=1:length(sub_regs)
  disp(j)
  sub_regs{j} = add_genes(D,sub_regs{j},cyto,urg,[1 0]);
end


%% Reconstruct focals data

s = size(D.dat,1);
n = size(D.dat,2);
%{
params = struct( ...
    'broad_or_focal','focal',...
    'broad_len_cutoff',broad_len_cutoff,...
    't_amp',thresh(1),'t_del',thresh(2),...
    'column_to_add',12,...
    'use_segarray',use_segarray,...
    'rows',size(D.dat,1),...
    'cols',size(D.dat,2) ...
);

verbose('reconstructing focal data',20);
focals = reconstruct_genomes(struct('amp',QA,'del',QD), params);


%% Write out list of peaks for each tumor type

verbose('writing amplification and deletion gene lists',20);

% loop over amplifications, then deletions
for k=1:2
  if k==1
    f = fopen([web_dir 'GCM_amp_regions_gene_lookup.' datext '.txt'],'w');
  else
    f = fopen([web_dir 'GCM_del_regions_gene_lookup.' datext '.txt'],'w');
  end
  fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Tissue type','Region','q-value',['Freq ' ...
                      'all'],'Freq focal','Freq high','Gene List','Summary Description');
  for j=1:length(sub_regs)
    disp([k j])
    cur_regs = sub_regs{j};
    cur_type = user_types{j};
    
    % Write summary string
    num_samples = sum(al_idx{j});
    sub_types = unique({D.sis(al_idx{j}).gcmtype});
    if length(sub_types) == 1
      summary_str = ['<p>This dataset (' char(cur_type) [') consists ' ...
                          'of '] num2str(num_samples) ' cancer samples.</p>'];
    else
      summary_str = ['<p>This dataset (' char(cur_type) ') consists of ' ...
                     num2str(num_samples)  ' cancer samples from ' ...
                     num2str(length(sub_types)) ' cancer subtypes: ' ...
                     '</p><br/><table><tr><th>Tissue Type</th><th># of Samples</th></tr> '];
      cur_idxs = cell(1,length(sub_types));
      for i=1:length(sub_types)
        cur_idxs{i} = al_idx{j}(strmatch(sub_types{i}, ...
                                     {D.sis(al_idx{j}).gcmtype},'exact'));
      end
      nsts = cellfun(@length,cur_idxs);
      [S si] = sort(nsts,'descend');
      
      for i=1:length(si)
        cur_idx = cur_idxs{si(i)};
        summary_str = [summary_str '<tr><td>' char(sub_types{si(i)}) '</td><td>' ...
                       num2str(length(cur_idx)) '</td></tr> '];
      end
      summary_str = [summary_str  '</table><br/>'];
    end
    ll = cellfun(@length,sub_regs{j});    
    summary_str = [summary_str '<p> ' num2str(ll(1)) ' peak regions of amplification and ', ...
                   num2str(ll(2)) ' peak regions of deletion were identified in this dataset.</p>'];
           
    for i=1:length(cur_regs{k})
      cur_reg = cur_regs{k}(i);
      q_value = cur_reg.resid_qv;
      peak_loc = ['chr' num2str(cur_reg.chrn) ':' ...
                  num2str(D.pos(cur_reg.peak_wide_st)) '-' ...
                  num2str(D.pos(cur_reg.peak_wide_en))];
      genes = cur_reg.genes;
      str = '';
      if ~isempty(genes)
        for l=1:length(genes)-1
          str = [str genes{l} '; '];
        end
        str = [str genes{end}];
      end
      
      xx1 = D.dat(cur_reg.peak,al_idx{j});
      if k==1
        xx2 = focals.amp(cur_reg.peak,al_idx{j});
      else
        xx1 = -xx1;
        xx2 = focals.del(cur_reg.peak,al_idx{j});
      end
      
      % calculate frequencies 
      freq_all = mean(xx1 >= thresh(k));
      freq_high = mean(xx1 >= high_thresh(k));
      freq_focal = mean(xx2 >= thresh(k));
 %{     
      % frequency all
      yy = xx1;
      yy(yy < 0.1) = 0;
      yy(yy > 0) = 1;
      freq_all = mean(yy);
      
      % frequency focal
      yy = xx2;
      yy(yy < 0.1) = 0;
      yy(yy > 0) = 1;
      freq_focal = mean(yy);
      
      % frequency high
      yy = xx1;
      yy(yy < 1) = 0;
      yy(yy > 0) = 1;
      freq_high = mean(yy);
%}      
      fprintf(f,'%s\t%s\t%1.3g\t%0.3f\t%0.3f\t%0.3f\t%s\t%s\n',cur_type,peak_loc,q_value,freq_all,freq_focal,freq_high,str,summary_str);
    end
  end
  fclose(f);
end % loop over amplifications,deletions

%% split peaks into amps and dels (amp_sub_regs and del_sub_regs) 
% construct lists of genes for each peak type (amp_genes and del_genes)
amp_sub_regs = cell(1,length(sub_regs));
del_sub_regs = cell(1,length(sub_regs));
amp_sub_q = cell(1,length(sub_regs));
del_sub_q = cell(1,length(sub_regs));
amp_genes = cell(1,length(sub_regs));
del_genes = cell(1,length(sub_regs));

% loop over disease types
for j=1:length(sub_regs)
  % split peaks into amps and dels
  cur_regs = sub_regs{j};
  amp_sub_regs{j} = cur_regs{1};
  del_sub_regs{j} = cur_regs{2};
  cur_q = sub_q{j};
  amp_sub_q{j} = cur_q{1};
  del_sub_q{j} = cur_q{2};
  if isempty(amp_sub_regs{j})
    xxa = {};
  else
    xxa = {amp_sub_regs{j}.genes};
  end
  if isempty(del_sub_regs{j})
    xxd = {};
  else
    xxd = {del_sub_regs{j}.genes};
  end
  % construct gene lists
  amp_genes{j} = {};
  del_genes{j} = {};
  for i=1:length(amp_sub_regs{j})
    amp_genes{j} = cat(2,amp_genes{j},xxa{i});
  end
  for i=1:length(del_sub_regs{j})
    del_genes{j} = cat(2,del_genes{j},xxd{i});
  end
end % loop over disease types


%% For each gene, make gene tables
%load([subtype_dir 'Breast/gene_stats.mat'])

% Create matrix representing relevant data for each gene in each type
% Matrix is (# types) x 9 x (# genes)
% 9 columns represent (range):
%  col 1: present (1) or absent from peak (0)
%  col 2: chrn # (1-23)
%  col 3: peak start (NaN or int)
%  col 4: peak end (NaN or int)
%  col 5: # genes in peak (0,Inf)
%  col 6: q-value (0,1)
%  col 7: frequency all (0,1)
%  col 8: frequency focal (0,1)
%  col 9: frequency high (0,1)

amp_data = NaN(ntypes,9,length(urg));
del_data = NaN(ntypes,9,length(urg));

verbose('creating individual gene tables for %d unique genes',20,length(urg));

%% loop over unique genes
for i=1:length(urg)
%!for i=1:length(ui)
  if mod(i,100) == 0
    disp(i)
  end
  gene = urg(i);
  amp_data(:,2,i) = gene.chrn;
  del_data(:,2,i) = gene.chrn;

  gene.midsnp = floor(median(gene.snps));

  % optimization: extract gene data before looping through types
  amp_snps = gene.snps;
  gene_maxdat = max(D.dat(amp_snps,:),[],1);
  gene_maxfoc = max(focals.amp(amp_snps,:),[],1);

  if ~isempty(gene.snps)
    %% loop over disease types
    for j=1:ntypes
      % check if gene is in an amp peak for this disease type
      if ~isempty(intersect(amp_genes{j},gene.symb))
        % gene in a peak
        amp_data(j,1,i) = 1;
        xx = cellfun(@(x) ~isempty(intersect(gene.symb,x)), ...
                     {amp_sub_regs{j}.genes},'UniformOutput',false);
        idx = find(arrayfun(@(l) xx{l} == 1,1:length(xx)));
        num_genes_per_idx = cellfun(@length,{amp_sub_regs{j}(idx).genes});
        % if more than one reg containing gene, pick smallest
        [m mi] = min(num_genes_per_idx);
        idx = idx(mi);
        amp_data(j,5,i) = m;
        cur_reg = amp_sub_regs{j}(idx);
        amp_data(j,3,i) = D.pos(cur_reg.peak_wide_st);
        amp_data(j,4,i) = D.pos(cur_reg.peak_wide_en);
      else
        % gene not in a peak
        amp_data(j,1,i) = 0;
        if isempty(amp_sub_regs{j})
          on_chr = [];
        else
          on_chr = find([amp_sub_regs{j}.chrn] == gene.chrn);
        end
        if ~isempty(on_chr)
          [m mi] = min(abs(D.pos([amp_sub_regs{j}(on_chr).peak])- ...
                           D.pos(gene.midsnp)));
          idx = on_chr(mi);
          amp_data(j,5,i) = length(amp_sub_regs{j}(idx).genes);
          cur_reg = amp_sub_regs{j}(idx);
          amp_data(j,3,i) = D.pos(cur_reg.peak_wide_st);
          amp_data(j,4,i) = D.pos(cur_reg.peak_wide_en);
        end
      end
      
      amp_data(j,6,i) = amp_sub_q{j}(gene.midsnp);
      
      % frequency all amplifications
      amp_data(j,7,i) = mean(gene_maxdat(al_idx{j}) >= thresh(1));
      % frequency focal amps
      amp_data(j,8,i) = mean(gene_maxfoc(al_idx{j}) >= thresh(1));
      % frequency high-level amps
      amp_data(j,9,i) = mean(gene_maxdat(al_idx{j}) >= high_thresh(1));
%{
      %Frequency all
      xx = full(max(D.dat(gene.snps,al_idx{j}),[],1));
      xx(xx < thresh) = 0;
      xx(xx > 0) = 1;
      amp_data(j,7,i) = mean(xx);
      
      %Frequency focal
      xx = full(max(focals.amp(gene.snps,al_idx{j}),[],1));
      xx(xx < thresh) = 0;
      xx(xx > 0) = 1;
      amp_data(j,8,i) = mean(xx);
      
      %Frequency high
      xx = full(max(D.dat(gene.snps,al_idx{j}),[],1));
      xx(xx < 1) = 0;
      xx(xx > 0) = 1;
      amp_data(j,9,i) = mean(xx);
 %}
    end
    
    %% deletion processing
    
    % precalculate widest slice of all gene isoforms for SegArray efficiency
    % (so different isoforms can be chosen for each disease type based on q-value)  
    del_slice = D.dat(gene.segs,:);
    focdel_slice = focals.del(gene.segs,:);

    
    % Find all isoforms of the gene that are in gene_gistic_rg
    isoforms = gene.iso_idx;
    
    % loop over disease types
    for j=1:ntypes
      % check if gene is in a del peak for this disease type
      if ~isempty(intersect(del_genes{j},gene.symb))
        % gene in a peak
        del_data(j,1,i) = 1;
        xx = cellfun(@(x) ~isempty(intersect(gene.symb,x)), ...
                     {del_sub_regs{j}.genes},'UniformOutput',false);
        idx = find(arrayfun(@(l) xx{l} == 1,1:length(xx)));
        num_genes_per_idx = cellfun(@length,{del_sub_regs{j}(idx).genes});
        % if more than one reg containing gene, pick smallest
        [m mi] = min(num_genes_per_idx);
        idx = idx(mi);
        del_data(j,5,i) = m;
        cur_reg = del_sub_regs{j}(idx);
        del_data(j,3,i) = D.pos(cur_reg.peak_wide_st);
        del_data(j,4,i) = D.pos(cur_reg.peak_wide_en);
      else
        % gene not in a peak
        del_data(j,1,i) = 0;
        if isempty(del_sub_regs{j})
          on_chr = [];
        else
          on_chr = find([del_sub_regs{j}.chrn] == gene.chrn);
        end
        if ~isempty(on_chr)
          [m mi] = min(abs(D.pos([del_sub_regs{j}(on_chr).peak])- ...
                           D.pos(gene.midsnp)));
          idx = on_chr(mi);
          del_data(j,5,i) = length(del_sub_regs{j}(idx).genes);
          cur_reg = del_sub_regs{j}(idx);
          del_data(j,3,i) = D.pos(cur_reg.peak_wide_st);
          del_data(j,4,i) = D.pos(cur_reg.peak_wide_en);
        end
      end
      
      % q-value: find isoform with minimum q-value for deletions
      if ~isempty(isoforms)
        [temp_qv min_isoform] = min([gene_regs{j}(isoforms).q]);
        iso_snps = gene.iso_segs{min_isoform};
%!      iso_snps = gene_regs{j}(isoforms(min_isoform)).snps - snp_base;
      else
        temp_qv = del_sub_q{j}(gene.midsnp);
        iso_snps = gene.snps - snp_base;
      end
      del_data(j,6,i) = temp_qv;
      
      % deletion and high deletion frequency
      mindels = min(del_slice(iso_snps,al_idx{j}),[],1);
      del_data(j,7,i) = mean(mindels <= -thresh(2));
      del_data(j,9,i) = mean(mindels <= -high_thresh(2));
      % focal deletion frequency
      minfocdels = max(focdel_slice(iso_snps,al_idx{j}),[],1);
      del_data(j,8,i) = mean(minfocdels <= -thresh(2));
    end
  end
end


%write intermediate files
save(save_analysis_file,'urg','user_types','keep','amp_data','del_data');

% Calculate how many peaks each gene is in
num_amp_sig = arrayfun(@(x) sum(amp_data(isubsets,6,x) <= .25),1: ...
                        size(amp_data,3));
num_del_sig = arrayfun(@(x) sum(del_data(isubsets,6,x) <= .25),1: ...
                        size(del_data,3));

num_amp_sig_peaks = arrayfun(@(x) sum(amp_data(isubsets,6,x) ...
                    <= .25 & amp_data(isubsets,1,x) == 1),1:size(amp_data,3));

num_del_sig_peaks = arrayfun(@(x) sum(del_data(isubsets,6,x) ...
                    <= .25 & del_data(isubsets,1,x) == 1),1:size(del_data,3));

hca = histc(num_amp_sig,0:ntypes);
hcd = histc(num_del_sig,0:ntypes);

%% loop over gene output files

for i=1:length(urg)
  if mod(i,250) == 0
    disp(i)
  end
  
  gene = urg(i);
  
  if ~isnan(amp_data(1,1,i))
    f = fopen([web_dir regexprep(gene.symb,'/','_') '.txt'],'w');
    fprintf(f,'%s\t%s %s%s%s%s%s%s%s\n',...
              'Result:',gene.symb,'(chr',deblank(num2str(gene.chrn)),...
              ':',deblank(num2str(gene.start)),'-',...
              deblank(num2str(gene.end)),')');
    
    % update gene summaries
    cur_mat = amp_data(:,:,i);
    fprintf(f,'\n%s:\t%s\n\n','Amplification Summary',...
              write_gene_summary(gene,size(D.dat,2),user_types,...
              cur_mat,num_amp_sig,num_amp_sig_peaks,1,isubsets,root_type_idx));
    
    % amplification summary
    fprintf(f,'\n\t\t\t\t\t%s\n\t\t\t\t\t\t\t\t%s\n','Amplifications','Frequency Amplified');
    fprintf(f,'\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Tumor subset','In a focal driver peak?','Nearest focal driver peak','Number of genes in peak','q-value','Overall','Focal','High-level');
    
    % sort order is all_cancers, then in-peak by q value, then out-peak by q-value
    in_peak = non_rootypes(find(cur_mat(non_rootypes,1)));
    out_peak = setdiff(non_rootypes,in_peak);
    [~,si1] = sort(cur_mat(in_peak,6),'ascend');
    [~,si2] = sort(cur_mat(out_peak,6),'ascend');
    si = [root_type_idx in_peak(si1) out_peak(si2)];
    
    dp_str = {'Yes','No'};
    for j=1:length(si)
      if ~isnan(cur_mat(si(j),3))
        peak_loc = ['chr' num2str(cur_mat(si(j),2)) ':' num2str(cur_mat(si(j),3)) ...
                    '-' num2str(cur_mat(si(j),4))];
      else
        peak_loc = 'No peak on chromosome';
      end
      fprintf(f,'\t%s\t\t%s\t%s\t%2.0f\t%1.3g\t%0.4f\t%0.4f\t%0.4f\n', ...
              char(user_types(si(j))),dp_str{2-cur_mat(si(j),1)},peak_loc,cur_mat(si(j),5:9));
    end
    
    % deletion summary
    cur_mat = del_data(:,:,i);
    fprintf(f,'\n%s:\t%s\n\n','Deletion Summary', ...
              write_gene_summary(gene,size(D.dat,2),user_types, ...
              cur_mat,num_del_sig,num_del_sig_peaks,0,isubsets,root_type_idx));
    fprintf(f,'\n\t\t\t\t\t%s\n\t\t\t\t\t\t\t\t%s\n','Deletions','Frequency Deleted');
    fprintf(f,'\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Tumor subset','In a focal driver peak?','Nearest focal driver peak','Number of genes in peak','q-value','Overall','Focal','High-level');
    
    % sort order is all_cancers, then in-peak by q value, then out-peak by q-value
    in_peak = find(cur_mat(non_rootypes,1))';
    out_peak = setdiff(non_rootypes,in_peak);
    [~,si1] = sort(cur_mat(in_peak,6),'ascend');
    [~,si2] = sort(cur_mat(out_peak,6),'ascend');
    si = [root_type_idx in_peak(si1) out_peak(si2)];
    
    dp_str = {'Yes','No'};
    for j=1:length(si)
      if ~isnan(cur_mat(si(j),3))
        peak_loc = ['chr' num2str(cur_mat(si(j),2)) ':' num2str(cur_mat(si(j),3)) ...
                    '-' num2str(cur_mat(si(j),4))];
      else
        peak_loc = 'No peak on chromosome';
      end
      fprintf(f,'\t%s\t\t%s\t%s\t%2.0f\t%1.3g\t%0.4f\t%0.4f\t%0.4f\n', ...
              char(user_types(si(j))),dp_str{2-cur_mat(si(j),1)},peak_loc,cur_mat(si(j),5:9));
    end
    
    fclose(f);
  end
end
%}

%% Print tumor summary file

% the file structure is a two rows, first row is sample, second row is
% tumor type

%!tumors = strmatch('EMPTY',celltypes);

Dt_gene = reduce_to_genes(D,urg,'symb',1);
%!Dt_gene = collapse_to_genes(D,urg,{urg.symb},'symb');

unique_ttypes = unique({Dt_gene.sis.gcmtype});

idx = [];
for j=1:length(unique_ttypes)
  idx = [idx; strmatch(unique_ttypes{j},{Dt_gene.sis.gcmtype},'exact')];
end

Dt_gene = reorder_D_cols(Dt_gene,idx);

Dt_gene.sis = struct('name',D.sdesc(idx)','gcmtype',{Dt_gene.sis.gcmtype});

%%

f = fopen([web_dir 'GCM_tumors_gene_summary.txt'],'w');

for i=1:size(Dt_gene.dat,2)
  fprintf(f,'\t%s',char(Dt_gene.sis(i).name));
end
fprintf(f,'\n');
for i=1:size(Dt_gene.dat,2)
  fprintf(f,'\t%s',char(Dt_gene.sis(i).gcmtype));
end
fprintf(f,'\n');
for j=1:size(Dt_gene.dat,1)
  modi(j,250)
  if exist([web_dir Dt_gene.gdesc{j} '.txt'],'file')
    fprintf(f,'%s',Dt_gene.gdesc{j});
    for i=1:size(Dt_gene.dat,2)
      fprintf(f,'\t%1.2f',2+Dt_gene.dat(j,i));
    end
    fprintf(f,'\n');
  end
end
  
fclose(f);

toc
verbose('Gene summary done',20);
