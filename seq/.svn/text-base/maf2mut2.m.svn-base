function maf2mut2(P)
% maf2mut2(P)
%
% reads maf files and writes mut file
%
% required parameters:
%   P.
%     project       ('TCGA' or 'TSP')
%     maf_dir
%     maf_fname
%     num_fields    (one for each maf)
%     mut_dir
%     mut_fname
%     patients_to_use_file
%
% Mike Lawrence 2008-04-09
%
% Reads MAF files, combines/vocab-normalizes/filters, and outputs MUT file.
%
% 1. Use MAF files as is, as long as there is a header row.
% 2. Automatically ignore "biotage" column from Baylor files
%    and "Sequencing_Phase" column from WU file.
% 3. Combine all files
% 4. Standardize all vocabulary sets
% 5. Sort rows by chromosome then start position
% 6. Check for redundant/conflicting records
% 7. Look up "phase" from all_targets.fixed_names.txt
% 8. Look up "strand" by matching to target
% 9. Filtering:  eliminate:
%    step1: patients marked "do not use" in master sample info file
%    step2: records with mutation_status of unknown/indeterminate/none
%    step3: sites that were germline/LOH in any sample, except valid+somatic
%    step4: sites seen in >1 sample (except 1 random rep), except valid+somatic
%    step5: sites with conflicting mutation_status or validation_status
%    --> USE DESCRIPTIVE NAMES FOR FILTERING REASONS, RATHER THAN NUMBERS
%
% Added 2008-04-16
%
% 1. Include "context" column listing reference sequence context +- 5bp.
% 2. Include "newbase" column listing what the mutation was a change to.
%
% Added 2008-04-29
%
% 1. Re-called mutations for synonymous/missense/nonsense
%    based on my frame calls in target list
%
% Added 2008-05-20
%
% 1. Modified to handle "Verification_status" column
%
% Added 2008-06-23
% 
% 1. Added "use" column for one-step filtering
%
% Added 2008-07-08
%
% 1. Added "hypermutated" column
% 2. Added reset of random number generator so result is the same on each run
%
% Modified 2008-10-16
%
% 1. Converted to callable function
%
% Modified 2008-10-29
%
% 1. Made filtering labels more informative for Germline/LOH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'project',[]);
P=impose_default_value(P,'maf_dir',[]);
P=impose_default_value(P,'maf_fname',[]);
P=impose_default_value(P,'num_fields',[]);
P=impose_default_value(P,'mut_dir',[]);
P=impose_default_value(P,'mut_fname',[]);
P=impose_default_value(P,'patients_to_use_file',[]);
P=impose_default_value(P,'hm_file','/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/hypermutated.txt');
P=impose_default_value(P,'targ_fname','/xchip/tcga/gbm/analysis/lawrence/cov/all_targets.fixed_names-WITH_FRAME.txt');
P=impose_default_value(P,'shr_fname','/xchip/tcga/gbm/analysis/lawrence/maf2mut/shared_genes.txt');
P=impose_default_value(P,'pause_after_each_step',false);
P=impose_default_value(P,'halt_on_first_nonstandard_vocabulary',false);
P=impose_default_value(P,'report_conflicts',false);
P=impose_default_value(P,'suppress_type_doublechecking',false);
P=impose_default_value(P,'output_mut_file',true);
P=impose_default_value(P,'include_original_columns',false);
P=impose_default_value(P,'include_dbSNP_info',false);
P=impose_default_value(P,'use_representative_of_repetitive',true);
P=impose_default_value(P,'use_validated_or_verified_mutations_only',true);
P=impose_default_value(P,'but_use_all_silent_mutations',true);
P=impose_default_value(P,'use_phase1_genes_only',false);
P=impose_default_value(P,'use_nonhypermutated_only',false);
P=impose_default_value(P,'random_seed',0);
P=impose_default_value(P,'check_class_assignments',false);

if isempty(P.project), error('Must supply P.project'); end
if isempty(P.maf_dir), error('Must supply P.maf_dir'); end
if isempty(P.maf_fname), error('Must supply P.maf_fname'); end
if isempty(P.num_fields), error('Must supply P.num_fields'); end
if isempty(P.mut_dir), error('Must supply P.mut_dir'); end
if isempty(P.mut_fname), error('Must supply P.mut_fname'); end
if isempty(P.patients_to_use_file), error('Must supply P.patients_to_use_file'); end

try

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LOAD FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Loading MAF files... ');

% mapping from structure fields to column names

if P.include_original_columns
% original columns
field = { ...
'Hugo_Symbol', 'Hugo_Symbol';...
'Entrez_Gene_Id', 'Entrez_Gene_Id';...
'Center', 'Center';...
'NCBI_Build', 'NCBI_Build';...
'Chromosome', 'Chrom';...
'Start_position', 'Start_position';...
'End_position', 'End_position';...
'Strand', 'Strand';...
'Variant_Classification', 'Variant_Classification';...
'Variant_Type', 'Variant_Type';...
'Reference_Allele', 'Reference_Allele';...
'Tumor_Seq_Allele1', 'Tumor_Seq_Allele1';...
'Tumor_Seq_Allele2', 'Tumor_Seq_Allele2';...
'dbSNP_RS', 'dbSNP_RS';...
'dbSNP_Val_Status', 'dbSNP_Val_Status';...
'Tumor_Sample_Barcode', 'Tumor_Sample';...
'Matched_Norm_Sample_Barcode', 'Match_Norm_Sample';...
'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele1';...
'Match_Norm_Seq_Allele2', 'Match_Norm_Seq_Allele2';...
'Tumor_Validation_Allele1', 'Tumor_Validation_Allele1';...
'Tumor_Validation_Allele2', 'Tumor_Validation_Allele2';...
'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele1';...
'Match_Norm_Validation_Allele2', 'Match_Norm_Validation_Allele2';...
'Verification_Status', 'Verification_Status';...
'Validation_Status', 'Validation_Status';...
'Mutation_Status', 'Mutation_Status';...
};
elseif P.include_dbSNP_info
field = { ...
'dbSNP_RS', 'dbSNP_RS';...
'dbSNP_Val_Status', 'dbSNP_Val_Status';...
};
else
  field = {};
end

field = [ field; { ...

% our columns
'patient', 'tumor_sample' ;...
'gene', 'hugo_symbol' ;...
'chr', 'chrom' ;...
'start', 'start_position' ;...
'end', 'end_position' ;...
'type', 'variant_classification' ;...
'status', 'mutation_status' ;...
'verification', 'verification_status' ;...
'validation', 'validation_status' ;...
'center', 'center' ;...

'ref', 'reference_allele' ;...
'ts1', 'tumor_seq_allele1' ;...
'ts2', 'tumor_seq_allele2' ;...
'mns1', 'Match_Norm_Seq_Allele1' ;...
'mns2', 'Match_Norm_Seq_Allele2' ;...
'tv1', 'Tumor_Validation_Allele1' ;...
'tv2', 'Tumor_Validation_Allele2' ;...
'mnv1', 'Match_Norm_Validation_Allele1' ;...
'mnv2', 'Match_Norm_Validation_Allele2' ;...
}];

mut=[];
mut.fno=[];
for fno=1:length(P.maf_fname)
   fname = fullfile(P.maf_dir,P.maf_fname{fno});
   if ~exist(fname,'file'), error('%s does not exist',fname); end
   maf = read_table(fname,repmat('%s',1,P.num_fields(fno)),char(9),1,'whitespace','\b\r');
    % extract the data we want
   for f=1:size(field,1)
     col = find(strncmpi(maf.headers{1},field{f,2},length(field{f,2})));
     if isempty(col), error('No %s column in %s', field{f,2}, P.maf_fname{fno}); end
     if fno==1   % start new field
        mut = setfield(mut, field{f,1}, maf.dat{col});
     else        % append
        mut = setfield(mut, field{f,1}, [getfield(mut, field{f,1}); maf.dat{col}]);
     end
   end
   mut.fno(end+1:end+length(maf.dat{1})) = fno;
end
mut.fno = mut.fno';

if P.include_dbSNP_info    % fill out blank fields
  for i=1:slength(mut)
    if isempty(mut.dbSNP_RS{i}), mut.dbSNP_RS{i}='-'; end
    if isempty(mut.dbSNP_Val_Status{i}), mut.dbSNP_Val_Status{i}='-'; end
  end
end

fprintf('done.\n');
if P.pause_after_each_step, keyboard; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ASSIGN MUTATIONS TO "CATEGORIES" BASED ON REF/NORM/TUM ALLELES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Assigning categories based on alleles... ');

nm = slength(mut);

mut.a = cell(nm,1);
for m=1:nm
  mut.a{m} = [mut.ref(m);mut.ts1(m);mut.ts2(m);mut.mns1(m);mut.mns2(m);...
             mut.tv1(m);mut.tv2(m);mut.mnv1(m);mut.mnv2(m)];
end

% replace "unknown"/"nd"/"af"/"N" with "?"
for m=1:nm
    for i=1:9
        if isempty(mut.a{m}{i}) || ...
        strcmpi(mut.a{m}{i},'unknown') || ...
        strcmpi(mut.a{m}{i},'nd') || ...
        strcmpi(mut.a{m}{i},'N') || ...
        strcmpi(mut.a{m}{i},'af')
             mut.a{m}{i} = '?';
        end
    end
end

% find point mutations
mut.ispm=true(nm,1);
for m=1:nm
    for i=1:9
        if mut.a{m}{i}(1) == '?', continue; end
        if mut.a{m}{i}(1) == '-', mut.ispm(m)=false; end
        if length(mut.a{m}{i})~=1, mut.ispm(m)=false; end
    end
end

% combine alleles into genotypes
% note: all are already uppercase
mut.ref_orig = mut.ref;
mut.ref=cell(nm,1);
mut.norm=cell(nm,1);
mut.tum=cell(nm,1);
mut.vnorm=cell(nm,1);
mut.vtum=cell(nm,1);
for m=1:nm
    mut.ref{m} = sort([mut.a{m}{1} mut.a{m}{1}]);    % reference genotype
    mut.norm{m} = sort([mut.a{m}{4} mut.a{m}{5}]);    % normal genotype
    mut.tum{m} = sort([mut.a{m}{2} mut.a{m}{3}]);    % tumor genotype
    mut.vnorm{m} = sort([mut.a{m}{8} mut.a{m}{9}]);    % normal genotype (validated)
    mut.vtum{m} = sort([mut.a{m}{6} mut.a{m}{7}]);    % tumor genotype (validated)
end

% now break into categories
%
%    indels
%
%    nonmutated:
%    ref AA = norm AA = tum AA
%
%    germline mutations:
%    ref AA -> norm AB = tum AB  heterozygous
%    ref AA -> norm BB = tum BB  homozygous
%    ref AA -> norm BC = tum BC  double-mutant
%    in general,
%    ref ~= (norm = tum)
%    
%    somatic mutations:
%    ref AA = norm AA -> tum AB  heterozygous
%    ref AA = norm AA -> tum BB  homozygous
%    ref AA = norm AA -> tum BC  double-mutant
%    in general,
%    (ref = norm) ~= tum
%    
%    somatic-upon-germline ("SUG") mutations
%    in general,
%    ref ~= norm ~= tum

mut.categ=cell(nm,1);
mut.newbase=repmat({'?'},nm,1);
for m=1:nm
  if ~mut.ispm(m)
     mut.categ{m} = 'Indel';
     mut.newbase{m} = 'Indel';
     continue;
  end
  ref = mut.ref{m};
  norm = mut.norm{m};
  tum = mut.tum{m};
  vnorm = mut.vnorm{m};
  vtum = mut.vtum{m};
  if ~strcmp(vtum,'??'), tum=vtum; end
  if ~strcmp(vnorm, '??'), norm=vnorm; end
  if all(ref==norm) && all(norm==tum)
    % unmutated
    mut.categ{m} = 'nonmut';
    mut.newbase{m} = ref(1);
  elseif ~all(ref==tum) && all(norm=='??')
     % unknown-status
     if any(ref==tum)
        mut.categ{m} = 'uhet';
        mut.newbase{m} = setdiff(tum,ref);
     elseif tum(1)==tum(2)
        mut.categ{m} = 'uhom';
        mut.newbase{m} = tum(1);
     else
        mut.categ{m} = 'udbl';
     end
  elseif ~all(ref==norm) && all(norm==tum)
     % germline  
     if any(ref==tum)
        mut.categ{m}='ghet';
        mut.newbase{m} = setdiff(tum,ref);
     elseif tum(1)==tum(2)
        mut.categ{m}='ghom';
        mut.newbase{m} = tum(1);
     else
        mut.categ{m}='gdbl';
     end
   elseif all(ref==norm) && ~all(norm==tum)
    % somatic
    if any(norm==tum)
       mut.categ{m}='shet';
       mut.newbase{m} = setdiff(tum,ref);
    elseif tum(1)==tum(2)
       mut.categ{m}='shom';
       mut.newbase{m} = tum(1);
    else
       mut.categ{m}='sdbl';
    end
 elseif ~all(ref==norm) && ~all(norm==tum)
    % somatic-upon-germline
    mut.categ{m}='SUG';
  else
    % other
    mut.categ{m}='other'
  end

end

fprintf('done.\n');
if P.pause_after_each_step, keyboard; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  STANDARDIZE VOCABULARY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Standardizing vocabulary... ');
bad_vocab_flag = false;

% allowed vocabulary for type

synonym = {
'Missense', 'Missense' ;...
'Missense', 'Missense_mutation' ;...
'Nonsense', 'Nonsense' ;...
'Nonsense', 'Non_Sense_Mutation' ;...
'Nonsense', 'Nonsense_Mutation' ;...
'Nonsense', 'Nonstop' ;...
'Nonsense', 'Nonstop_Mutation' ;...
'Splice_site', 'Splice_site_SNP' ;...
'Splice_site', 'Spice_Site_SNP' ;...
'Splice_site', 'Split_Site_SNP' ;...
'Splice_site', 'Splice_site' ;...
'Splice_site', 'Splice_region' ;...
'Splice_site', 'Splice site' ;...
'Splice_site', 'Exon_boundary' ;...
'Splice_site', 'Splice_Region_Del' ;...
'Splice_site', 'Splice_Region_Ins' ;...
'Splice_site', 'Silent_Splice' ;...
'Splice_site', 'Splice_Site_Mutation' ;...
'Inframe_Ins', 'Inframe_Ins' ;...
'Frameshift_Del', 'Frame_Shift_Del' ;...
'Splice_site', 'Splice_Site_Del' ;...
'Inframe_Del', 'Inframe_Del' ;...
'Indel', 'Indel' ;...
'Frameshift','frameshift' ;...
'Frameshift_Del', 'Frame_Shift_ Del' ;...
'Inframe_Ins', 'In_Frame_Ins' ;...
'Frameshift', 'frame_shift' ;...
'Frameshift_Ins', 'Frame_Shift_Ins' ;...
'Frameshift_Ins', 'Frameshift_ins' ;...
'Inframe_Del', 'in_frame_deletion' ;...
'Inframe_Ins', 'in_frame_insertion' ;...
'Frameshift_Del', 'Frameshift_del' ;...
'Splice_site', 'Splice_Site_Indel' ;...
'Inframe_Del', 'In_Frame_ Del' ;...
'Inframe_Del', 'In_Frame_Del' ;...
'Synonymous', 'Synonymous' ;...
'Synonymous', 'Silent' ;...
'Synonymous', 'Silent_Mutation' ;...
'Targeted_Region', 'Targeted_region' ;...
'Targeted_Region', '3_prime_flanking_region' ;...
'Targeted_Region', '3_prime_untranslated_region' ;...
'Targeted_Region', '5_prime_flanking_region' ;...
'Targeted_Region', '5_prime_untranslated_region' ;...
'Targeted_Region', 'Intronic' ;...
'Unknown', 'Unknown' ;...
'Unknown', 'undefined' ;...
'Unknown', 'Consensus_Error' ;...
'Unknown', 'None' ;...
'Unknown', '' ;...
};

for i=1:slength(mut)
    v = find(strcmpi({synonym{:,2}}',mut.type{i}));
    if isempty(v)
      if P.halt_on_first_nonstandard_vocabulary
        error('Invalid type %s in %s', mut.type{i}, P.maf_fname{mut.fno(i)});
      else
        fprintf('Invalid type %s in %s\n', mut.type{i}, P.maf_fname{mut.fno(i)});
        bad_vocab_flag = true;
      end
    else
      mut.type{i} = synonym{v,1};
    end
end

% allowed vocabulary for status

synonym = {
'Germline', 'germ line' ;...
'Germline', 'Germline' ;...
'LOH', 'LOH' ;...
'LOH', '"somatic,loh"' ;...
'Indeterminate', 'Indeterminate' ;...
'Unknown', 'Unknown' ;...
'Unknown', '' ;...
'Unknown', 'na' ;...
'Unknown', 'None' ;...
'Unknown', 'Variant' ;...
'Unknown', 'SomaticAlleleMismatch' ;...
'Somatic', 'Somatic' ;...
'Somatic', '"somatic, heterozygous"' ;...
'Somatic', '"somatic, homozygous"' ;...
'Somatic', '"somatic,heterozygous"' ;...
'Somatic', '"somatic,homozygous"' ;...
'Valid', 'Valid';...
};

for i=1:slength(mut)
    v = find(strcmpi({synonym{:,2}}',mut.status{i}));
    if isempty(v)
      if P.halt_on_first_nonstandard_vocabulary
        error('Invalid status %s in %s', mut.status{i}, P.maf_fname{mut.fno(i)});
      else
        fprintf('Invalid status %s in %s\n', mut.status{i}, P.maf_fname{mut.fno(i)});
        bad_vocab_flag = true;
      end
    else
      mut.status{i} = synonym{v,1};
    end
end

% allowed vocabulary for validation

synonym = {
'Unknown', 'Unknown' ;...
'Valid', 'Valid' ;...
'Wildtype', 'Wildtype' ;...
'Valid', 'Biotage' ;...
'Unknown', '-' ;...
'Unknown', '' ;...
'Unknown', 'Disregard manual review' ;...
'Unknown', 'allelemismatch' ;...
};

for i=1:slength(mut)
    v = find(strcmpi({synonym{:,2}}',mut.validation{i}));
    if isempty(v)
      if P.halt_on_first_nonstandard_vocabulary
        error('Invalid validation %s in %s', mut.validation{i}, P.maf_fname{mut.fno(i)});
      else
        fprintf('Invalid validation %s in %s\n', mut.validation{i}, P.maf_fname{mut.fno(i)});
        bad_vocab_flag = true;
      end
    else
      mut.validation{i} = synonym{v,1};
    end
end

% allowed vocabulary for verification

synonym = {
'Germline', 'Germline' ;...
'Somatic', 'Somatic' ;...
'Unknown', 'Unknown' ;...
'Valid', 'Valid' ;...
'Wildtype', 'Wildtype' ;...
'Unknown', '-' ;...
'Unknown', '' ;...
'Unknown', 'U' ;...
'Unknown', 'Disregard manual review' ;...
'Valid', 'Verified' ;...
};

for i=1:slength(mut)
    v = find(strcmpi({synonym{:,2}}',mut.verification{i}));
    if isempty(v)
      if P.halt_on_first_nonstandard_vocabulary
        error('Invalid verification %s in %s', mut.verification{i}, P.maf_fname{mut.fno(i)});
      else
        fprintf('Invalid verification %s in %s\n', mut.verification{i}, P.maf_fname{mut.fno(i)});
        bad_vocab_flag = true;
      end
    else
      mut.verification{i} = synonym{v,1};
    end
end

if bad_vocab_flag
  fprintf('Quitting due to nonstandard vocabulary');
  keyboard
end

fprintf('done.\n');
if P.pause_after_each_step, keyboard; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  OTHER SIMPLE CONVERSIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% truncate patient IDs to first 12 characters

if strcmp(P.project,'TCGA')
for i=1:slength(mut)
    mut.patient{i} = mut.patient{i}(1:12);
end
end

% construct "coords" column = chrX:start-end

mut.coords = cell(slength(mut),1);
for i=1:slength(mut)
    mut.coords{i} = ['chr' mut.chr{i} ':' mut.start{i} '-' mut.end{i}];
end

% construct "site" column = gene|chrX:stard-end

mut.site = cell(slength(mut),1);
for i=1:slength(mut)
   mut.site{i} = [mut.gene{i} '|' mut.coords{i}];
end

% convert "start" and "end" from strings to numbers

mut.start = str2double(mut.start);
mut.end = str2double(mut.end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CHECK FOR REDUNDANT / CONFLICTING RECORDS
%     look for records with same patient + coords
%     for matching records,
%       compare type, status, validation, and center
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Scanning for redundant/conflicting records... ');

[p pi pj] = unique(mut.patient);
[c ci cj] = unique(mut.coords);

discard_idx = [];
for i=1:slength(mut)
    if ~mod(i,1000), fprintf('%d/%d ',i,slength(mut)); end
    if find(discard_idx==i), continue; end    % skip records that have already been discarded
    x = find(pj==pj(i));
    x = intersect(x, find(cj==cj(i)));
    if length(x)>1                      % found records with same patient+coords

        % mark the other records as to-be-discarded
        discard_idx = [discard_idx; setdiff(x,i)];

        % concatenate center names
        tmp = '';
        u = unique(mut.center(x));
        for j=1:length(u)
            if length(tmp), tmp = [tmp ', ']; end
            tmp = [tmp u{j}];
        end
        mut.center{i} = tmp;

        % check for conflicts in newbase

        u = unique(mut.newbase(x));
        if length(u)>1
          u = setdiff(u, '?');
          if length(u)==1   % only one is non-?
            mut.newbase{i} = u{1}; 
          else
            mut.newbase{i} = '?';
            if P.report_conflicts
              fprintf('Newbase conflict in %s %s %s\n', mut.patient{i}, mut.gene{i}, mut.coords{i});
              [mut.center(x) mut.newbase(x)]
            end
          end
        end

        % check for conflicts in type, status, validation            

        if length(unique(mut.type(x)))>1
            if P.report_conflicts && length(unique(mut.center(x)))==1
              fprintf('Type conflict in %s %s %s\n', mut.patient{i}, mut.gene{i}, mut.coords{i});
%              [mut.center(x) mut.ref(x) repmat({'->'},length(x),1) mut.tum1(x) ...
%               repmat({'+'},length(x),1) mut.tum2(x) mut.type(x)]
            end
            tmp = '';
            u = unique(mut.type(x));
            for j=1:length(u)
               if length(tmp), tmp = [tmp '/']; end
               tmp = [tmp u{j}];
            end
            mut.type{i} = tmp;
        end

        if length(unique(mut.status(x)))>1
            if P.report_conflicts
              fprintf('Status conflict in %s %s %s\n', mut.patient{i}, mut.gene{i}, mut.coords{i});
              [mut.center(x) mut.status(x)]
            end
            tmp = '';
            u = unique(mut.status(x));
            for j=1:length(u)
               if strcmp(u{j}, 'Unknown'), continue; end    % ignore "Unknowns"
               if length(tmp), tmp = [tmp '/']; end
               tmp = [tmp u{j}];
            end
            mut.status{i} = tmp;
        end

        if length(unique(mut.validation(x)))>1
           if P.report_conflicts && length(unique(mut.center(x)))==1
              fprintf('Validation conflict in %s %s %s\n', mut.patient{i}, mut.gene{i}, mut.coords{i});
              [mut.center(x) mut.validation(x)]
           end
           tmp = '';
            u = unique(mut.validation(x));
            for j=1:length(u)
               if strcmp(u{j}, 'Unknown'), continue; end    % ignore "Unknowns"
                if length(tmp), tmp = [tmp '/']; end
               tmp = [tmp u{j}];
            end
            mut.validation{i} = tmp;
        end

        if length(unique(mut.verification(x)))>1
           if P.report_conflicts && length(unique(mut.center(x)))==1
              fprintf('Verification conflict in %s %s %s\n', mut.patient{i}, mut.gene{i}, mut.coords{i});
              [mut.center(x) mut.verification(x)]
           end
           tmp = '';
            u = unique(mut.verification(x));
            for j=1:length(u)
               if strcmp(u{j}, 'Unknown'), continue; end    % ignore "Unknowns"
                if length(tmp), tmp = [tmp '/']; end
               tmp = [tmp u{j}];
            end
            mut.verification{i} = tmp;
        end
    end
end

% remove discarded records

keep_idx = setdiff([1:slength(mut)], discard_idx);
mut = reorder_struct(mut,keep_idx);

fprintf('done.\n');
fprintf('Discarded %d redundant records; kept %d/%d.\n',...
  length(discard_idx),length(keep_idx),length(discard_idx)+length(keep_idx));
if P.pause_after_each_step, keyboard; end

nmut = slength(mut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SORT RECORDS: by chromosome, then by start position
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Sorting records... ');

mut.chrtmp = regexprep(mut.chr, '^(\d)$', '0$1');      % prepend single-digit chromosomes with 0 for sort
mut.sorttmp = cell(nmut,1);
for i=1:nmut
    mut.sorttmp{i} = sprintf('%s:%09d', mut.chrtmp{i}, mut.start(i));
end

mut = sort_struct(mut,'sorttmp');

fprintf('done.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MARK SHARED GENES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(P.project,'TCGA')
  shared_genes = read_table(P.shr_fname, '%s%f', char(9), 1, 'whitespace', '\b\r'); 
  mut.shared = repmat({'No'},nmut,1);
  for j=1:length(shared_genes.dat{1})
    x = find(strcmp(shared_genes.dat{1}{j}, mut.gene));
    mut.shared(x) = repmat({'Yes'}, length(x), 1);
  end
else
  mut.shared = repmat({'?'}, nmut, 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MATCH RECORDS TO TARGET LIST; EXTRACT PHASE AND STRAND
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(P.project, 'TCGA')

fprintf('Loading target list... ');
target = read_table(P.targ_fname, '%s%s%s%f%f%s%s%f', char(9), 1, 'Whitespace', '\b\r');
target.dat{1} = regexprep(target.dat{1},'\..*','');    % remove spliceform identifiers (eg .2)
fprintf('done.\n');

fprintf('Matching mutations to targets... ');

mut.phase = cell(nmut,1);
mut.strand = cell(nmut,1);
mut.targets = cell(nmut,1);

[g gi gj] = unique(target.dat{1});
for i=1:nmut
    gidx = find(strcmp(mut.gene{i},g));
    gtargs = [];
    targs = [];
    if ~isempty(gidx)
      gtargs = find(gj==gidx);
      ts = target.dat{4}(gtargs);
      te = target.dat{5}(gtargs);
      targs = find(mut.start(i)<te & mut.end(i)>ts);   % find which target it's in
      if ~isempty(targs)
        mut.targets{i} = gtargs(targs);
      end
    end
    if ~isempty(gtargs)         % if matched to a gene
       u = unique(target.dat{7}(gtargs));  % phase
       if length(u)>1, error('Conflict over phase!'); end
       mut.phase{i} = u{1};
       u = unique(target.dat{6}(gtargs));  % strand
       if length(u)>1
         u = unique(target.dat{6}(gtargs(targs)));
         u = setdiff(u,{'.'});
         if isempty(u)
           u = {'?'};
         elseif length(u)>1
           fprintf('Conflict over strand in gene %s\n',g{gidx}); u, keyboard
         end
       end
       if isempty(u), fprintf('Problem with strand calling\n'); keyboard, end
       mut.strand{i} = u{1};
    else                        % failed to match to any gene:  "unknown", "LOC115648"
       mut.phase{i} = 'Unknown';
       mut.shared{i} = 'Unknown';
       mut.strand{i} = 'Unknown';
    end
end

fprintf('done.\n');
if P.pause_after_each_step, keyboard; end

else  % TSP
  mut.strand = repmat({'?'}, nmut, 1);
  mut.phase = repmat({'?'}, nmut, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  TAKE COMPLEMENT OF NEWBASE FOR (-)STRAND MUTATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nmut
   if strcmp(mut.strand{i}, '-')
      if ismember(mut.newbase{i}, 'ACGT')
        mut.newbase{i} = my_seqrcomplement(mut.newbase{i});
      end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  IDENTIFY CONTEXT AND MUTATION CLASS
%
%     first mark indels (len>1)
%     then for non-indels, look up genome region +- 5 bp and
%         1. assign "class" (i.e. context)
%     and 2. re-call "type" (synonymous/missense/nonsense)
%            (necessary because some centers did it wrong)
%            For mutations in genes with multiple spliceforms,
%            break down "type" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Identifying context and mutation class... ');

mut.class = cell(nmut,1);
mut.context = cell(nmut,1);
mut.recalled = false(nmut,1);
mut.oldtype = mut.type;

indels = grep('In|Del|Frameshift', mut.type, 1);
mut.class(indels) = repmat({'Indel'},length(indels),1);
mut.context(indels) = repmat({'Indel'},length(indels),1);
pointmuts = setdiff([1:nmut],indels);

for p=1:length(pointmuts)    % for each point mutant

   if ~mod(p,2000), fprintf('%d/%d ', p, length(pointmuts)); end
   i = pointmuts(p);

   % Look up context of mutation

   context = lower(genome_region(mut.chr{i}, mut.start(i)-5, mut.start(i)+5));
   if mut.strand{i} == '-', context = my_seqrcomplement(context); end
   context(6) = upper(context(6));
   mut.context{i} = context;

   if strcmp(P.project, 'TCGA')

   % Manually exclude TCGA-02-0010|CDKN2A|chr9:21961053-21961053 from double-checking
   
   if strcmp(mut.patient{i},'TCGA-02-0010') && strcmp(mut.gene{i},'CDKN2A') && (mut.start(i)==21961053)
     immune_to_doublecheck = true;
   else
     immune_to_doublecheck = false;
   end

   % Double-check mutation type (syn/mis/non)
   %   if newbase is known and mutation was matched to a target

   if ~strcmp(mut.newbase(i),'?') && ~strncmp('Splice',mut.type{i},6) && ~isempty(mut.targets{i})
     oldtype = mut.type{i};
     ntarg = length(mut.targets{i});
     t=1;   % only consider first target for now
     targ = mut.targets{i}(t);
     frame = target.dat{8}(targ);
     if frame~=1000   % if frame is known for this target
       pos = mut.start(i);
       ts = target.dat{4}(targ);
       te = target.dat{5}(targ);
       strand = target.dat{6}{targ};
        newbase = mut.newbase{i};
       if strcmp(strand,'-')
         bp_into_target = te - pos;
       else
         bp_into_target = pos - ts;
       end
       mutframe = mod(frame + bp_into_target, 3);
       newcontext = context;
       if length(newbase)~=1, newbase = '?'; end
       newcontext(6:6) = newbase;
       if mutframe==0, codon=context(6:8); newcodon=newcontext(6:8);
       elseif mutframe==1, codon=context(5:7); newcodon=newcontext(5:7);
       else codon=context(4:6); newcodon=newcontext(4:6); end
       if strcmpi(codon, newcodon)
         newtype = 'Non-mutation';
       else
         oldaa = my_nt2aa(codon);
         newaa = my_nt2aa(newcodon);
         if newaa==oldaa
           newtype = 'Synonymous';
         elseif newaa=='*'
           newtype = 'Nonsense';
         else
           newtype = 'Missense';
         end
       end
       % if Type needs to be recalled on the basis of the doublecheck:
       if ~P.suppress_type_doublechecking
         if ~strcmpi(oldtype,newtype)
            if immune_to_doublecheck
              fprintf('%s|%s|chr%s:%d-%d immune to type recall (would be %s->%s)\n',...
                mut.patient{i},mut.gene{i},mut.chr{i},mut.start(i),mut.end(i),...
                mut.type{i},newtype);
            else
               mut.recalled(i) = true;
               mut.type{i} = newtype;
            end
         end
       end
     end
   end  % end of mutation-type doublechecking

   else   %  TSP
   end

   % FINALLY: assign mutation-context class

   dna = upper(context(5:7));
   if dna(2) == 'A' || dna(2) == 'T'
       mut.class{i} = dna(2);
   elseif dna(2) == 'C'
       if dna(3) == 'G'
            mut.class{i} = 'C (CpG)';
       elseif dna(1) == 'T'
            mut.class{i} = 'C (TpC)';
       else
            mut.class{i} = 'C (other)';
       end
   elseif dna(2) == 'G'
       if dna(1) == 'C'
            mut.class{i} = 'G (CpG)';
       elseif dna(3) == 'A'
            mut.class{i} = 'G (GpA)';
       else
            mut.class{i} = 'G (other)';
       end
   else
            mut.class{i} = 'Unknown';
   end
end

mut.typechange = cell(nmut,1);
for i=1:nmut
   mut.typechange{i} = concat({mut.oldtype{i},mut.type{i}},' -> ');
end

fprintf('done.\n');

if P.suppress_type_doublechecking
  fprintf('Type doublechecking was suppressed.\n');
else
  fprintf('%d mutations were re-called for Type.\n', sum(mut.recalled));
  count(mut.typechange(mut.recalled));
end

if P.pause_after_each_step, keyboard; end

if P.check_class_assignments
  fprintf('Checking class assignments... ');
  for i=1:nmut
      if strcmp('Indel',mut.class{i}), continue; end
      r = mut.ref{i};
      if mut.strand{i} == '-', r = my_seqrcomplement(r); end
      if r ~= mut.class{i}(1)
          fprintf('Class conflict: %d %s %s %s %s\n', i, mut.type{i}, mut.strand{i}, r, mut.class{i});
      end
  end
end
fprintf('done.\n');
if P.pause_after_each_step, keyboard; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ASSIGN "Supertype"
%
%  Supertype            included Types
%  --------------       -----------------------------------------
%  Non-mutation         Non-mutation
%  Synonymous           Synonymous
%  Missense             Missense
%  Nonsense             Nonsense
%  Unknown              Unknown
%  Other                Indel, Splice_site, Targeted_region, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mut.supertype = repmat({'Other'},nmut,1);

idx = grep('^(Non-mutation|Synonymous|Missense|Nonsense|Unknown)$', mut.type, 1);
mut.supertype(idx) = mut.type(idx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ASSIGN "CALL" and "EVIDENCE"
%
%  Based on Validation_status, Verification_status, and Mutation_status
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mut.call = cell(nmut,1);
mut.evidence = cell(nmut,1);
for m=1:nmut
  ver = mut.verification{m};
  val = mut.validation{m};
  stat = mut.status{m};
  if strcmp(val, 'Valid')
    evid = 'Validated';
    call = stat;
  elseif strcmp(val, 'Wildtype')
    evid = 'Validated';
    call = val;
  elseif strcmp(ver, 'Valid')
    evid = 'Verified';
    call = stat;
  elseif grep('Germline|Somatic|Wildtype',ver,1)
    evid = 'Verified';
    call = ver;
  else
    evid = 'Sequenced';
    call = stat;
  end
  mut.evidence{m} = evid;
  mut.call{m} = call;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FILTERING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Filtering records... ');

mut.filtered = repmat({'OK'},nmut,1);

if strcmp(P.project, 'TCGA')

% FILTERING STEP 0:  Remove non-mutations and Call=Wildtype

nonmuts = find(strcmp(mut.supertype, 'Non-mutation'));
wildtype = find(strcmp(mut.call, 'Wildtype'));
nonmuts = union(nonmuts, wildtype);
mut.filtered(nonmuts) = repmat({'Non-mutation'}, length(nonmuts), 1);

% FILTERING STEP 1:  Remove patients that do not appear in the list of patients to use

tmp = load_struct(P.patients_to_use_file);
patients_to_use = tmp.name;

for i=1:nmut
    x = find(strcmp(patients_to_use, mut.patient{i}));
    if isempty(x)
       mut.filtered{i} = 'Do_not_use';
    end
end

% FILTERING STEP 2:  Remove records with Call=Unknown

for i=1:nmut
    if ~strcmp(mut.filtered{i}, 'OK'), continue; end   % skip already-filtered
    if strcmp(mut.call{i}, 'Unknown')
       mut.filtered{i} = 'Status_unknown';
    end
end

% FILTERING STEP 3:  Remove sites that are call=germline/LOH in any sample,
%                    but always spare somatic+(valid/verified)

nonfiltered = find(strcmp(mut.filtered,'OK'));

% find which records are germline or LOH

germline = grep('Germline', mut.call, 1);
LOH = grep('LOH',mut.call, 1);

% find which records are Evidence=(Verified/Validated) + Call=Somatic

vv = grep('Validated|Verified', mut.evidence, 1);
somatic = find(strcmp(mut.call,'Somatic'));
vv_somatic = intersect(vv,somatic);

% for each unique site, count how many samples have germline or LOH

[nms,nmsi,nmsj] = unique(mut.site);
for s=1:length(nms)
    a = find(nmsj==s);
    a = intersect(a, nonfiltered);   % only consider nonfiltered records
    % records listed in "a" have this site
    bg = intersect(a, germline);
    bg = intersect(bg, nonfiltered);
    bl = intersect(a, LOH);
    bl = intersect(bl, nonfiltered);
                                     % only consider nonfiltered records
    % records listed in "b" have this site and are germline or LOH
    % if there are any records in "b", then we need to filter out all
    % records in "a" except those with "VALID/VERIFIED SOMATIC"
    if length(bg)>0 && length(bl)>0
        c = setdiff(a, vv_somatic);
        for i=1:length(c)
            if ~isempty(grep('Germline|LOH',mut.call{c(i)},1))
              mut.filtered{c(i)} = mut.call{c(i)};
            else
              mut.filtered{c(i)} = 'Germline/LOH_in_other_sample';
            end
        end
    elseif length(bg)>0
        c = setdiff(a, vv_somatic);
        for i=1:length(c)
            if ~isempty(grep('Germline|LOH',mut.call{c(i)},1))
              mut.filtered{c(i)} = mut.call{c(i)};
            else
              mut.filtered{c(i)} = 'Germline_in_other_sample';
            end
        end
    elseif length(bl)>0
        c = setdiff(a, vv_somatic);
        for i=1:length(c)
            if ~isempty(grep('Germline|LOH',mut.call{c(i)},1))
              mut.filtered{c(i)} = mut.call{c(i)};
            else
              mut.filtered{c(i)} = 'LOH_in_other_sample';
            end
        end
    end

end

% FILTERING STEP 4:  Remove sites seen in >1 sample (except one randomly chosen representative)
%                    but always spare (verified/valid)+somatic

nonfiltered = find(strcmp(mut.filtered,'OK'));     % update nonfiltered list

rand('state',P.random_seed);   % so the result is deterministic

for s=1:length(nms)
    a = find(nmsj==s);
    a = intersect(a, nonfiltered);  % only consider nonfiltered records
    % records listed in "a" have this site
    if (length(a)>1)
        % if more than one record in "a", then:
        %     if NONE are valid somatic,
        %             choose one at random and keep it, while filtering out
        %             all others.
        %     else (at least one is valid/verified somatic)
        %             keep all the valid/verified somatic ones, while filtering out
        %             all others.
        b = intersect(a, vv_somatic);
        if (isempty(b))  % no records in "a" are valid_somatic
            keeper = ceil(rand(1) * length(a));
            for i=1:length(a)
                if (i==keeper)
                    mut.filtered{a(i)} = 'Representative_of_repetitive';
                else
                    mut.filtered{a(i)} = 'Repetitive';
                end
            end
        else   % at least one is valid somatic
            c = setdiff(a, vv_somatic);
            for i=1:length(c)
                mut.filtered{c(i)} = 'Repetitive';
            end
        end
    end
end

% FILTERING STEP 5:  Records with conflicts in Call or Validation

for i=1:nmut
    if ~strcmp(mut.filtered{i}, 'OK'), continue; end   % skip already-filtered
    if ~isempty(find(mut.call{i}=='/'))
       mut.filtered{i} = 'Call_conflict';
    elseif ~isempty(find(mut.validation{i}=='/'))
       mut.filtered{i} = 'Validation_conflict';
    end
end


else    % TSP
end


fprintf('done.\n');
if P.pause_after_each_step, keyboard; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FINALLY, add "use" field with "Do_not_use", "Use_silent", or "Use_nonsilent"
% and "hypermutated" field with "hypermutated" or "non-hypermutated"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hm = load_struct(P.hm_file);

if strcmp(P.project,'TCGA')
  x=grep('TCGA-06-0150',mut.patient,1);
  mut.filtered(x) = repmat({'Do_not_use'},length(x),1);
  x=find(strcmp('CHEK2|chr22:27413910-27413910',mut.site) & ~strcmp('Synonymous',mut.type));
  mut.evidence(x) = repmat({'Sequenced'},length(x),1);
end

mut.use = repmat({'Do_not_use'},nmut,1);

% apply usage criteria

if P.use_representative_of_repetitive
  use = grep('^OK|Representative_of_repetitive$',mut.filtered,1);
else
  use = grep('^OK$',mut.filtered,1);
end

if P.use_phase1_genes_only
  p1 = grep('1',mut.phase,1);
  use = intersect(use,p1);
end

if isfield(hm,'name'), tmp=hm.name;
elseif isfield(hm,'patient'),tmp=hm.patient;
else tmp={}; end
h = find(ismember(mut.patient,tmp));

if P.use_nonhypermutated_only
  use = setdiff(use,h);
end

vv = grep('^Verified|Validated$',mut.evidence,1);
usevv = intersect(use,vv);

sil = grep('Synonymous',mut.type,1);

if P.use_validated_or_verified_mutations_only && ~P.but_use_all_silent_mutations
  use_silent = intersect(usevv,sil);
else
  use_silent = intersect(use,sil);
end

if P.use_validated_or_verified_mutations_only
  use_nonsilent = setdiff(usevv,sil);
else
  use_nonsilent = setdiff(use,sil);
end

mut.use(use_silent) = repmat({'Use_silent'},length(use_silent),1);
mut.use(use_nonsilent) = repmat({'Use_nonsilent'},length(use_nonsilent),1);

mut.hypermutated = repmat({'Non-hypermutated'},nmut,1);
mut.hypermutated(h) = repmat({'Hypermutated'},length(h),1);

if P.pause_after_each_step, keyboard; end

if P.output_mut_file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  OUTPUT MUT FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mapping from output columns to field names

if P.include_original_columns
outfield = { ...

% original columns
'Hugo_Symbol', 'Hugo_Symbol';...
'Entrez_Gene_Id', 'Entrez_Gene_Id';...
'Center', 'Center';...
'NCBI_Build', 'NCBI_Build';...
'Chromosome', 'Chromosome';...
'Start_position', 'Start_position';...
'End_position', 'End_position';...
'Strand', 'Strand';...
'Variant_Classification', 'Variant_Classification';...
'Variant_Type', 'Variant_Type';...
'Reference_Allele', 'Reference_Allele';...
'Tumor_Seq_Allele1', 'Tumor_Seq_Allele1';...
'Tumor_Seq_Allele2', 'Tumor_Seq_Allele2';...
'dbSNP_RS', 'dbSNP_RS';...
'dbSNP_Val_Status', 'dbSNP_Val_Status';...
'Tumor_Sample_Barcode', 'Tumor_Sample_Barcode';...
'Matched_Norm_Sample_Barcode', 'Matched_Norm_Sample_Barcode';...
'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele1';...
'Match_Norm_Seq_Allele2', 'Match_Norm_Seq_Allele2';...
'Tumor_Validation_Allele1', 'Tumor_Validation_Allele1';...
'Tumor_Validation_Allele2', 'Tumor_Validation_Allele2';...
'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele1';...
'Match_Norm_Validation_Allele2', 'Match_Norm_Validation_Allele2';...
'Verification_Status', 'Verification_Status';...
'Validation_Status', 'Validation_Status';...
'Mutation_Status', 'Mutation_Status';...
};
elseif P.include_dbSNP_info
outfield = { ...
'dbSNP_RS', 'dbSNP_RS';...
'dbSNP_Val_Status', 'dbSNP_Val_Status';...
};
else
  outfield = {};
end

outfield = [outfield; { ...
% our columns

%%%%% NOTE: FIRST FIVE FIELDS MUST BE KEPT AS FOLLOWS FOR IGV:

'chr', 'chr' ;...
'start', 'start' ;...
'end', 'end' ;...
'patient', 'patient' ;...
'type', 'type' ;...

%%%%% OTHER FIELDS ARE ALLOWED TO VARY.

'gene', 'gene' ;...
'phase', 'phase' ;...
'center', 'center' ;...
'shared', 'shared' ;...
'strand', 'strand' ;...
'site', 'site' ;...
'supertype', 'supertype' ;...

'ref', 'ref_orig' ;...
'ts1', 'ts1' ;...
'ts2', 'ts2' ;...
'mns1', 'mns1' ;...
'mns2', 'mns2' ;...
'tv1', 'tv1' ;...
'tv2', 'tv2' ;...
'mnv1', 'mnv1' ;...
'mnv2', 'mnv2' ;...

'context', 'context' ;...
'class', 'class' ;...
'newbase', 'newbase' ;...
'status', 'status' ;...
'validation', 'validation' ;...
'call', 'call' ;...
'evidence', 'evidence' ;...
'filtered', 'filtered' ;...
'use', 'use' ;...
'hypermutated', 'hypermutated' ;...
}];

fprintf('Outputting MUT file... ');

M=[];
for i=1:length(outfield)
   M = setfield(M,outfield{i,1},getfield(mut,outfield{i,2}));
end

save_struct(M,[P.mut_dir P.mut_fname]);

fprintf('done.\n');

if P.pause_after_each_step, keyboard; end

end % if P.output_mut_file


catch me, excuse(me); end
