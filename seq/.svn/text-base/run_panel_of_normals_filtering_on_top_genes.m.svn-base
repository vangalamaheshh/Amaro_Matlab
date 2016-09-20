function run_panel_of_normals_filtering_on_top_genes(maf,genelist,bamlist,outdir,P)
% run_panel_of_normals_filtering_on_top_genes(maf,genelist,bamlist,outdir,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','');
P = impose_default_value(P,'dynamically_update_bamlist',false);

if isempty(P.build)
  fprintf('Assuming hg19.\n');
  P.build = 'hg19';
end

if nargin<4, error('need four parameters'); end
if ~ischar(maf) || ~ischar(genelist) || ~ischar(bamlist) || ~ischar(outdir)
  error('four parameters should be strings');
end

demand_files({maf;genelist;bamlist});
ede(outdir);

fprintf('Loading bamlist\n');
bam = load_lines(bamlist);
if contains(bam{1},char(9))
  bam = load_struct(bamlist);
  if isfield(bam,'name') && ~isfield(bam,'bam'), bam = rename_field(bam,'name','bam'); end
  if isfield(bam,'bamfile') && ~isfield(bam,'bam'), bam = rename_field(bam,'bamfile','bam'); end
else
  if strcmp(bam{1},'name') || strcmp(bam{1},'bam') || strcmp(bam{1},'bamfile')
    bam = bam(2:end);
  end
  tmp=bam; 
  bam=[]; bam.bam=tmp; clear tmp
end

fprintf('Loading genelist\n');
gene = load_lines(genelist);
if contains(gene{1},char(9))
  gene = load_struct(genelist);
  if isfield(gene,'gene') && ~isfield(gene,'name'), gene = rename_field(gene,'gene','name'); end
else
  if strcmp(gene{1},'name') || strcmp(gene{1},'gene')
    gene = gene(2:end);
  end
  tmp=gene;
  gene=[]; gene.name=tmp; clear tmp
end

fprintf('Loading mutations\n');
m = load_struct(maf);
m = add_simple_fieldnames(m);

if P.dynamically_update_bamlist
  fprintf('Updating bamlist\n');
  bam = update_bamlist(bam);
end

istop = ismember(m.gene,gene.name);
top = reorder_struct(m,istop);
btm = reorder_struct(m,~istop);
fprintf('Will screen %d/%d mutations.\n',slength(top),slength(m));

save_struct(top,[outdir '/topgenes.maf']);
save_struct(bam,[outdir '/bams_to_screen.txt']);

survey_panel_of_normals_for_mutations([outdir '/topgenes.maf'],bam.bam,outdir,P);

topnew = load_struct([outdir '/topgenes.PoN_filtered.maf']);
mnew = concat_structs_keep_all_fields({topnew,btm});

[path name ext] = fileparts(maf);
newmaf = fullfile(outdir,[name '.PoN_filtered' ext]);
save_struct(mnew,newmaf);
fprintf('\nSaved final filtered MAF to %s\n',newmaf);


