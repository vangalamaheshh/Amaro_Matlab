function M = load_all_mutation_data3(P,isetname)

if iscell(P)
  % multi-set load mode
  if exist('isetname','var'), error('use P.isetname for multi-set load mode'); end
  M = cell(length(P),1);
  for i=1:length(P), M{i} = load_all_mutation_data3(P{i}); end
  return
end

if ~exist('isetname','var'), isetname = []; end

if ~exist('P','var'),P=[];end
P=impose_default_value(P,'targlist','*required*');
P=impose_default_value(P,'summed_cov_track','*required*');
P=impose_default_value(P,'build','*required*');
P=impose_default_value(P,'isetname',isetname);

demand_file(P.targlist);
demand_file(P.summed_cov_track);

P.ignore_coverage = true;
M = load_all_mutation_data2(P,isetname);

M.cov = [];
M.cov.targ =load_target_file(P.targlist);
M.cov.targ.gidx = listmap(M.cov.targ.gene,M.gene.name);

ng = slength(M.gene);
for g=1:ng
  idx = find(M.cov.targ.gidx==g);
  M.gene.chr(g,1) = M.cov.targ.chr(idx(1));
  M.gene.start(g,1) = min(M.cov.targ.start(idx));
  M.gene.end(g,1) = max(M.cov.targ.end(idx));
  M.gene.len(g,1) = sum(M.cov.targ.len(idx));
  if isfield(M.cov.targ,'gc'), M.gene.gc(g,1) = weighted_mean(M.cov.targ.gc(idx),M.cov.targ.len(idx)); end
  if isfield(M.cov.targ,'type'), M.gene.type(g,1) = M.cov.targ.type(idx(1)); end
end

M.file.summed_cov_track = P.summed_cov_track;
M.build = P.build;


