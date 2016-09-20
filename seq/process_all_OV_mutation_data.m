function process_all_OV_mutation_data(P)
% Mike Lawrence 2009-07-02

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'ABI_patient_list_input_file','*required*');
P=impose_default_value(P,'capture_patient_list_input_file','*required*');
P=impose_default_value(P,'ABI_mutation_list_input_file','*required*');
P=impose_default_value(P,'capture_mutation_list_input_file','*required*');
P=impose_default_value(P,'indel_mutation_list_input_file','*required*');
P=impose_default_value(P,'combined_mutation_list_output_file','*required*');
P=impose_default_value(P,'combined_patient_list_output_file','*required*');

try

% LOAD PATIENT LISTS

fprintf('Loading patient lists\n');
ABI_patients = load_struct(P.ABI_patient_list_input_file);
CAP_patients = load_struct(P.capture_patient_list_input_file);

% LOAD ABI MUTATIONS

fprintf('Loading ABI mutations\n');
ABI = load_struct(P.ABI_mutation_list_input_file);
ABI.dataset = repmat({'ABI'},slength(ABI),1);
ABI.chr = str2double(ABI.chr);
ABI = rename_field(ABI,{'ref','ts1','ts2'},{'ref_allele','tum_allele1','tum_allele2'});

% LOAD CAPTURE MUTATIONS

fprintf('Loading capture mutations\n');
CAP = load_struct(P.capture_mutation_list_input_file);
% remove any excess header lines resulting from the file concatenation
CAP = reorder_struct(CAP,setdiff(1:slength(CAP),find(strcmp(CAP.build,'build'))));
CAP = reorder_struct(CAP,grep('Missense|Synonymous|Nonsense|Splice|Read-through',CAP.type,1));
CAP.chr = convert_chr(CAP.chr);
for i=1:slength(CAP),CAP.patient{i,1} = CAP.tumor_barcode{i}(1:12); end
CAP.dataset = repmat({'CAP'},slength(CAP),1);

% LOAD INDELS

fprintf('Loading indels\n');
INDEL = load_indel_data(P);

% COMBINE ALL MUTATION DATA

fprintf('Combining all mutation data\n');
flds = {'patient','chr','start','end','ref_allele','tum_allele1','tum_allele2','type','gene','dataset'};
M1 = concat_structs({keep_fields(ABI,flds),keep_fields(CAP,flds),keep_fields(INDEL,flds)});

% unique by patient+site
for i=1:slength(M1), id{i,1} = [M1.patient{i} '|chr' M1.chr(i) ':' M1.start{i} '-' M1.end{i}]; end
[u ui uj] = unique(id);
M2 = reorder_struct(M1,ui);
for i=1:length(u),M2.dataset{i} = concat(M1.dataset(find(uj==i)),'+');end

% ANNOTATE

fprintf('Annotating mutations\n');
idx = grep('Missense|Nonsense|Read-through|Synonymous',M2.type,1);
M3 = reorder_struct(M2,idx); M3.pos = str2double(M3.start);
M3 = annotate_sites(M3);

M = M2;
nm = slength(M);
M.transcript = repmat({'-'},nm,1);
M.proteinchange = repmat({'-'},nm,1);
for i=1:length(idx)
  if ismember(M3.type{i},{'Missense','Nonsense','Read-through','Synonymous'})
    M.transcript{idx(i)} = M3.transcript{i};
    M.proteinchange{idx(i)} = M3.proteinchange{i};
  else
    M.transcript{idx(i)} = '?';
    M.proteinchange{idx(i)} = '?';
  end
end

% ADD OTHER COLUMNS

fprintf('Adding additional information\n');
M.filtered = repmat({'OK'},slength(M),1);
M.type = regexprep(M.type,'Read-through','Nonsense');   % (temporary fix)

% mutation classes (categories)
M.class = upper(M.ref_allele);
M = make_numeric(M,{'start','end'});
idx = union(grep('In|Del',M.type,1),find(M.end-M.start>0));
M.class(idx) = repmat({'Indel'},length(idx),1);

% splice-site indels --> frameshift     (temporary fix)
idx = find(strcmp(M.type,'Splice_site') & strcmp(M.class,'Indel'));
M.type(idx) = repmat({'Frameshift'},length(idx),1);

% annotate "class" (base categories)
for i=1:slength(M)
  if strcmp(M.class{i},'C')
    base = upper(genome_region(M.chr(i),M.start(i)+1));
    if strcmp(base,'G'), M.class{i} = 'C (CpG)';
    else M.class{i} = 'C (other)'; end
  elseif strcmp(M.class{i},'G')
    base = upper(genome_region(M.chr(i),M.start(i)-1));
    if strcmp(base,'C'), M.class{i} = 'G (CpG)';
    else M.class{i} = 'G (other)'; end
  end
end

% change DNPs to single mutations

fprintf('Identifying DNPs\n');
n = nan(slength(M),1);
idx = cell(slength(M),1);
pointmuts = grep('Synonymous|Nonsense|Missense',M.type,1);
for i=1:slength(M)
  idx{i} = intersect(pointmuts,find(strcmp(M.patient,M.patient(i)) & ...
     M.chr==M.chr(i) & abs(M.start - M.start(i))<=3));
  n(i) = length(idx{i});
end
keep = 1:slength(M);
M4 = M;

for i=1:slength(M4), j=idx{i};
  if n(i)>1
    if ismember('Nonsense',M.type(j)), M4.type{i} = 'Nonsense';
    else M4.type{i} = 'Missense'; end
    M4.start(i) = min(M.start(j));
    M4.end(i) = max(M.end(j));
    ref = genome_region(M.chr(i),M4.start(i),M4.end(i));
    tum1_is_ref = strcmp(M.tum_allele1{i},M.ref_allele{i});
    tum2_is_ref = strcmp(M.tum_allele2{i},M.ref_allele{i});
    M4.ref_allele{i} = ref;
    tum = ref;
    for k=1:length(j); q = j(k);
      if ~strcmp(M.tum_allele1{q},M.ref_allele{q}), newseq = M.tum_allele1{q};
      else newseq = M.tum_allele2{q}; end
      tum([M.start(q):M.end(q)]-M4.start(i)+1) = newseq;
    end
    if ~tum1_is_ref, M4.tum_allele1{i} = tum; else M4.tum_allele1{i} = ref; end
    if ~tum2_is_ref, M4.tum_allele2{i} = tum; else M4.tum_allele2{i} = ref; end
    M4.proteinchange{i} = '?';
    keep = setdiff(keep, j(2:end));
  end
end
M = reorder_struct(M4,keep);

% add "site"
chr2 = convert_chr_back(M.chr);
for i=1:slength(M)
  M.site{i,1} = [M.gene{i} '|' chr2{i} ':' ...
    num2str(M.start(i)) '-' num2str(M.end(i))];
end


% SAVE COMBINED MUT FILE

fprintf('Saving combined mutation file %s\n',P.combined_mutation_list_output_file);
% first five columns in canonical order
X = keep_fields(M,{'chr','start','end','patient','type',...
  'class','site','gene','ref_allele','tum_allele1','tum_allele2',...
  'transcript','proteinchange','dataset','filtered'});
X.chr = convert_chr_back(X.chr);

save_struct(X,P.combined_mutation_list_output_file);


% SAVE COMBINED PATIENT FILE

fprintf('Saving combined patient file %s\n',P.combined_patient_list_output_file);
PAT = concat_structs({ABI_patients,CAP_patients});
[tmp ui uj] = unique(PAT.short);
PAT = reorder_struct(PAT,ui);
save_struct(PAT,P.combined_patient_list_output_file);

catch me; excuse(me); end
