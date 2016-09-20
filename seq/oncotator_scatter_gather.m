function m = oncotator_scatter_gather(infile,outfile,P);

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','hg19');
P = impose_default_value(P,'queue','hour');
P = impose_default_value(P,'mutations_per_job',3000);

if ~strcmp(P.build,'hg19'), error('only hg19 supported'); end

if exist(outfile,'file'), fprintf('Output file already exists\n'); return; end

% load and validate input file
fprintf('Loading input file.\n');
m = load_struct(infile);
flds = {'build','chr','start','end','ref_allele','tum_allele1','tum_allele2','tumor_barcode','normal_barcode'};
demand_fields(m,flds);
keep_fields(m,flds);
nm = slength(m);
njobs = ceil(nm/P.mutations_per_job);
if any(grepm('23|24',m.chr)), error('must specify X or Y, not 23 or 24'); end
if any(grepm('36|18',m.build)), error('only hg19 supported'); end

% divide input file into chunks
fprintf('Dividing into per-job files: ');
tempdir = [outfile '.tempdir']; ede(tempdir);
firstm = 1;
part_in_file = cell(njobs,1);
part_out_file = cell(njobs,1);
for i=1:njobs,if~mod(i,100), fprintf('%d/%d ',i,njobs); end
  lastm = min(nm,firstm+P.mutations_per_job-1);
  mi = reorder_struct(m,[firstm:lastm]);
  part_in_file{i,1} = [tempdir '/input.part' num2str(i) '.maf'];
  part_out_file{i,1} = [tempdir '/output.part' num2str(i) '.maf'];
  if ~exist(part_in_file{i},'file'), save_struct(mi,part_in_file{i}); end
  firstm = lastm+1;
end, fprintf('\n');

results = cell(njobs,1);
done = false;
while(~done)
  cmds={}; banners={};
  jobs = [];
  for i=1:njobs
    if ~exist(part_in_file{i},'file'), error('what?'); end
    if exist(part_out_file{i},'file'), continue; end
    cmd = ['source /broad/software/scripts/useuse;unuse Python-2.6;reuse .python-2.7.1-sqlite3-rtrees;reuse .zlib-1.2.6;'...
           '/xchip/tcga/Tools/oncotator/onco_env/bin/oncotator '...
           '-v --override_config=/xchip/cga1/annotation/db/tcgaMAFManualOverrides.config '...
           part_in_file{i} ' ' part_out_file{i} ' hg19 --no-multicore'];
    cmd = ['"' cmd '"'];
    banner = ['Onco' num2str(i)];
    cmds{end+1} = cmd; banners{end+1} = banner;
    if length(cmds)>=60, jobs = [jobs; bsub(cmds,banners,P)]; cmds={}; banners={}; end
  end
  if ~isempty(cmds), jobs = [jobs; bsub(cmds,banners,P)]; cmds={}; banners={}; end
  if ~isempty(jobs), bwait(jobs); end
  fprintf('All Oncotator jobs finished\n');

  % make sure all have the correct number of lines
  fprintf('Loading output files: \n');
  ok = false(njobs,1);
  for i=1:njobs, if ~mod(i,10), fprintf('%d/%d ',i,njobs); end
    if ok(i), continue; end % already loaded in previous pass
    if ~exist(part_out_file{i},'file'), fprintf('Output file %d not found!',i); ok(i) = false; continue; end
    results{i} = load_struct(part_out_file{i});
    expected_length = P.mutations_per_job;
    if i==njobs, expected_length = mod(nm,P.mutations_per_job); end
    if slength(results{i})~=expected_length, fprintf('File %d not correct length.\n',i); ok(i) = false; continue; end
    ok(i) = true;
  end

  if all(ok)
    done = true;
  else
    bad = find(~ok);
    fprintf('%d/%d files are bad/missing\n',length(bad),njobs);
    % easiest solution is to just delete them and start over
    fprintf('Deleting the truncated files... ');
    for i=1:length(bad), if exist(part_out_file{bad(i)}), delete(part_out_file{bad(i)}); end, end
    fprintf('and resubmitting the jobs\n');
  end
end

results = concat_structs(results);
fprintf('Writing final output file\n');
save_struct(results,outfile);






