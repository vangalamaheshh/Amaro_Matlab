function X = dRanger_blast_filter_process_results(jobid,X,outdir,tempdir,max_queries_per_blast_batch,blast_radius,max_hits)

if ~exist('jobid','var'), error('jobid is required'); end
if ~exist('X','var'), error('X is required'); end
if ~exist('outdir','var'), error('outdir is required'); end

if ~exist('max_queries_per_blast_batch','var'), max_queries_per_blast_batch = 25; end
if ~exist('blast_radius','var'), blast_radius = 250; end
if ~exist('max_hits','var'), max_hits = 10; end
if ~exist('tempdir','var'), tempdir = '/xchip/tcga_scratch/lawrence/tmp'; end

fprintf('dRanger_blast_filter_process_results\n\toutdir = %s\n\tjobid = %s\n',outdir,jobid);

nx = slength(X);

fprintf('Preprocessing blast results\n');
fnames = [tempdir '/blast_' jobid '_hit_batch_*'];
allfile = [tempdir '/blast_' jobid '_hits_all'];

%keyboard   %%%% <----------------------

system(['cat ' fnames ' | grep -vP "random|hap|^#" | sed -e ''s/^jump//'' -e ''s/\.end/\t/'' -e ''s/chr//'' '...
   '-e ''s/X/23/'' -e ''s/Y/24/'' -e ''s/M/0/'' > ' allfile]);

fprintf('Loading results\n');
H = load_struct(allfile,repmat('%f',1,13),0);
H = rename_field(H,{'col1','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11','col12','col13'},...
                 {'jump','end','chr','pctid','matchlen','n_mm','n_gaps','qstart','qend','hstart','hend','E', ...
                  'score'});
H.qstrand = repmat({'+'},slength(H),1);
H.hstrand = H.qstrand;

% replace start>end with explicit strand indicators

idx = find(H.qstart>H.qend);
if ~isempty(idx)
  tmp = H.qstart(idx); H.qstart(idx) = H.qend(idx); H.qend(idx) = tmp;
  H.qstrand(idx) = repmat({'-'},length(idx),1);
end

idx = find(H.hstart>H.hend);
if ~isempty(idx)
  tmp = H.hstart(idx); H.hstart(idx) = H.hend(idx); H.hend(idx) = tmp;
  H.hstrand(idx) = repmat({'-'},length(idx),1);
end

% analyze results

X.hits1 = nan(nx,1);
X.hits2 = nan(nx,1);
X.filterB = nan(nx,1);

fprintf('Analyzing results\n');
for i=1:nx
    fprintf('Jump %d --> ', i);
    chr1 = X.chr1(i); pos1 = round((X.min1(i)+X.max1(i))/2);
    chr2 = X.chr2(i); pos2 = round((X.min2(i)+X.max2(i))/2);
    mn1 = pos1-blast_radius; mx1 = pos1+blast_radius;
    mn2 = pos2-blast_radius; mx2 = pos2+blast_radius;
    jidx = find(H.jump==i);
    H1 = reorder_struct(H,jidx(H.end(jidx)==1));
    H2 = reorder_struct(H,jidx(H.end(jidx)==2));
    fprintf('(%d+%d) ',slength(H1),slength(H2));
    X.hits1(i) = slength(H1); X.hits2(i) = slength(H2);
    if slength(H1)==0 | slength(H2)==0
      fprintf('BLAST failed\n'); X.filterB(i) = NaN; continue;
    end
    idx1 = find(H1.chr==chr2 & H1.hstart<mx2 & H1.hend>mn2);
    idx2 = find(H2.chr==chr1 & H2.hstart<mx1 & H2.hend>mn1);
    if ~isempty(idx1) | ~isempty(idx2)
      fprintf('Rejected because of homology\n'); X.filterB(i) = 2; continue;
    end
    if slength(H1)>max_hits | slength(H2)>max_hits
      fprintf('Rejected because of too many hits\n'); X.filterB(i) = 1; continue;
    end
    fprintf('Could not be rejected\n'); X.filterB(i) = 0;
end

save_struct(X,[outdir '/dRanger_blast_filtered_' jobid '_X.txt']);
