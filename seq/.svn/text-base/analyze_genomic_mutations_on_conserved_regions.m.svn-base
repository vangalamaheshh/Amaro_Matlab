function analyze_genomic_mutations_on_conserved_regions(regfile,mutfile,covfile,mutrate,outstem)
% analyze_genomic_mutations_on_conserved_regions(regfile,mutfile,covfile,mutrate,outstem)
%
% if covfile is [], will use Nterr as Ncov

demand_file({regfile;mutfile});
if ~isempty(covfile), demand_file(covfile); end

fprintf('Loading mutations\n');
X = load_struct(mutfile);
X = require_fields_with_convert(X,{'indiv','chr','start'},{'patient','','',''});
if ~isnumeric(X.chr), X.chr=convert_chr(X.chr); end
X = make_numeric(X,{'start'});
X = sort_struct(X,{'chr','start'});

fprintf('Loading regions\n');
tmp = load(regfile); f = fieldnames(tmp); R = getfield(tmp,f{1});
require_fields(R,{'chr','start','end'});
if ~isfield(R,'gene'), R.gene = repmat({'?'},slength(R),1); end
if ~isfield(R,'zone'), R.zone = repmat({'?'},slength(R),1); end
if ~isfield(R,'ig'), R.ig = nan(slength(R),1); end
if ~isnumeric(R.chr), R.chr=convert_chr(R.chr); end
R = make_numeric(R,{'start','end','ig'});
if ~issorted(R.chr), fprintf('Sorting R\n'); R = sort_struct(R,{'chr','start','end'}); end

if ~isnumeric(mutrate) || mutrate<1e-8 || mutrate>1e-3, error('inappropriate mutrate'); end

try
  save_textfile('test',[outstem '.top1000_regions.txt']);
catch me
  error('can''t write output using outstem="%s"',outstem);
end

% remove gene annotations from non-Refseq genes
idx = find(strcmp(R.gene,'?'));
R.gene(idx) = repmat({'-'},length(idx),1);
R.zone(idx) = repmat({'IGR'},length(idx),1);

% mark "bad" genome areas (this includes the hypermutated Ig regions)
fprintf('Marking "bad" genome areas\n');
R.goodbad = get_context(R.chr,round((R.end+R.start)/2),'/xchip/cga1/lawrence/db/goodbad','.goodbad.mat');

fprintf('Extracting coverage\n');
% get coverage in each region
if ~isempty(covfile)
  R = extract_region_coverage_from_wig(R,covfile);
else
  fprintf('Assuming full coverage ');
  R.len = R.end-R.start+1;
  nsamps = length(unique(X.indiv));
  fprintf('with %d patients\n',nsamps);
  R.cov = R.len*nsamps;
end

fprintf('Counting mutations\n');
% count mutations in each region
R = match_mutations_to_regions(R,X);    % assumes X and R are both sorted

% remove mutations in uncovered regions
idx=find(R.cov==0 & R.nmuts>0);
R.nmuts(idx) = 0; R.nsamps(idx) = 0; R.details(idx) = repmat({''},length(idx),1);

fprintf('Computing significance\n');
% calculate significance
R.p = 1-binocdf(R.nmuts-1,R.cov,mutrate);
R.q = calc_fdr_value(R.p);
R = sort_struct(R,{'q','p','nsamps'},[1 1 -1]);

fprintf('Writing results\n');
save_struct(reorder_struct(R,1:1000),[outstem '.top1000_regions.txt']);
save([outstem '.all_regions.mat'],'R');

% remove Ig regions
fprintf('Writing results with Ig loci masked\n');
R = reorder_struct(R,~R.ig);
R.q = calc_fdr_value(R.p);
R = sort_struct(R,{'q','p','nsamps'},[1 1 -1]);

save_struct(reorder_struct(R,1:1000),[outstem '.top1000_regions.Igmasked.txt']);
save([outstem '.all_regions.Igmasked.mat'],'R');
