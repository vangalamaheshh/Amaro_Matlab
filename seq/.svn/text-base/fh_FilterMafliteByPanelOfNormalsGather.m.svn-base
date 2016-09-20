function fh_FilterMafliteByPanelOfNormalsGather(maf, bamlist, min_num_alt_reads, min_alt_to_ref_ratio, min_num_normals, sample_id)

fprintf('fh_FilterMafliteByPanelOfNormalsGather...\n');
demand_file(maf);
demand_file(bamlist);

if ~isnumeric(min_num_alt_reads), min_num_alt_reads=str2double(min_num_alt_reads); end
if ~isnumeric(min_alt_to_ref_ratio), min_alt_to_ref_ratio=str2double(min_alt_to_ref_ratio); end
if ~isnumeric(min_num_normals), min_num_normals=str2double(min_num_normals); end

% move files out of scatter job directories into survey_files directory
survey_dir = './survey_files'; ensure_dir_exists(survey_dir);
status = system(['mv scatter*/survey_file_*.txt ' survey_dir]);
if status~=0
  error('file move failed');
end

% see if maflite is empty
d = dir(maf);
if d.bytes==0
  fprintf('maf is empty\n');
  save_textfile('',[sample_id '.filtering_details.txt']);
  save([sample_id '.filtering_details.mat']);
  save_textfile('',[sample_id '.maf']);
  return
end


% load maflite
M = load_struct(maf);
orig_flds2 = {'build','chrom','start_position','end_position','ref_allele','alt_allele1','alt_allele2','tumor_name','normal_name'};
orig_flds = {'build','contig','start_position','end_position','ref_allele','alt_allele1','alt_allele2','tumor_name','normal_name'};
flds = {'build','chr','start','end','ref_allele','tum_allele1','tum_allele2','tumor_barcode','normal_barcode'};
for i=1:length(flds)
  if isfield(M,orig_flds{i}) && ~isfield(M,flds{i})
    M = rename_fields(M,orig_flds{i},flds{i});
  end
  if isfield(M,orig_flds2{i}) && ~isfield(M,flds{i})
    M = rename_fields(M,orig_flds2{i},flds{i});
  end
end

M = order_fields_first(M,flds);
M_input = M;
M = make_numeric(M,'start');
M.chr = convert_chr(M.chr);
M.newbase = find_newbase(M);

% load results from scatter jobs
B = load_lines(bamlist);
outfiles = regexprep(num2cellstr(1:length(B)),'(.*)',[survey_dir '/survey_file_$1.txt']);
fprintf('Loading results: ');
nn = 0;
ncols = get_colcount(outfiles{1},0);
if ncols~=7 && ncols~=9, error('wrong column count in first result file'); end
fmt = ['%f%f%s' repmat('%f',1,ncols-3)];
C = nan(slength(M),ncols-3,length(B));
for f=1:length(B), fprintf('%d/%d ',f,length(B));
  if exist(outfiles{f})
    tmp = read_table(outfiles{f},fmt,char(9),0);
    z = cat(2,tmp.dat{4:end});
    if size(z,1)~=slength(M), error('wrong number of mutations in results file %s',outfiles{f}); end
    C(:,:,f) = z;
    nn = nn + 1;
  else
    if ~exist(B{f})
      fprintf('Warning: BAM file not found: %s\n',B{f});
    else
      fprintf('Warning: No results for BAM file: %s\n',B{f});
    end
  end
end, fprintf('\n');

% count ref/alt alleles   (so far, ignoring del+ins)
bases = 'ACGT';
D = nan(slength(M),4,length(B));
for i=1:slength(M)
  mutref = M.ref_allele{i};
  mutalt = M.newbase{i};
  if length(mutref)~=1 || length(mutalt)~=1, continue; end
  ref = find(bases == M.ref_allele{i});
  alt = find(bases == M.newbase{i});
  if isempty(ref) || isempty(alt), continue; end
  others = setdiff(1:4,[ref alt]);
  oth1 = others(1);
  oth2 = others(2);
  D(i,:,:) = C(i,[ref alt oth1 oth2],:);
end

% count evidence
num_alt_reads = squeeze(D(:,2,:));
if slength(M)==1, num_alt_reads = num_alt_reads'; end
num_ref_reads = squeeze(D(:,1,:));
if slength(M)==1, num_ref_reads = num_ref_reads'; end
num_tot_reads = squeeze(sum(D,2));
if slength(M)==1, num_tot_reads = num_tot_reads'; end
alt_to_ref_ratio = num_alt_reads./num_ref_reads;
has_evidence = (num_alt_reads>=min_num_alt_reads & alt_to_ref_ratio>=min_alt_to_ref_ratio);
num_samps_with_evidence = sum(has_evidence,2);
keep = find(num_samps_with_evidence<min_num_normals);

% record details of evidence
M.normal_panel_evidence_details = repmat({'---'},slength(M),1);
B2 = regexprep(B,'.*/([^/]+)$','$1');
if length(unique(B2)) < length(unique(B)), B2=B; end
B3 = regexprep(B2,'(.*)-Normal.bam','$1');
if length(unique(B3)) < length(unique(B)), B3=B2; end
for i=1:slength(M)
  x = '';
  for j=1:size(B,1)
    if has_evidence(i,j)
      n = [B3{j} '(' num2str(num_alt_reads(i,j)) '/' num2str(num_tot_reads(i,j)) ')'];
      if isempty(x), x=n; else x = [x ';' n]; end
  end,end
  if ~isempty(x), M.normal_panel_evidence_details{i} = x; end
end
M.num_normals_in_panel = nn*ones(slength(M),1);
M.num_normals_with_evidence_of_alt_allele = num_samps_with_evidence;
M.remove_based_on_normals_panel = (num_samps_with_evidence>=min_num_normals);

% filter and save
save_struct(M,[sample_id '.filtering_details.txt']);
save([sample_id '.filtering_details.mat']);

M_filtered = reorder_struct(M_input,keep);
save_struct(M_filtered,[sample_id '.maf']);

fprintf('Done.\n');
